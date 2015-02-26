# kdArrange v.003
import bmesh
import bpy
import mathutils
from mathutils import Vector
from collections import defaultdict

from BioBlender.table_values import (values_fi, scale_vdw)


## helper functions

def bmesh_from_pyverts(verts):
    bm = bmesh.new()
    add_vert = bm.verts.new
    bm_verts = [add_vert(co) for co in verts]
    bm.verts.index_update()
    return bm


def build_ktree(v):
    # documentation/blender_python_api_2_70_release/mathutils.kdtree.html
    size = len(v)
    kd = mathutils.kdtree.KDTree(size)

    for i, vtx in enumerate(v):
        kd.insert(Vector(vtx), i)
    kd.balance()
    return kd


def get_vcol_layer(obj):
    vcols = obj.data.vertex_colors
    if not ('fi_cols' in vcols):
        vcol_layer = obj.data.vertex_colors.new('fi_cols')
    else:
        vcol_layer = vcols.get('fi_cols')
    return vcol_layer


def get_specific_dict(filepath):
    '''
    This generates a dict from the .pdb to remap the specific elements in a chain
    to it's generic element type.
    {'OD2': 'O', 'CD1': 'C', 'CH2': 'C', 'CG': 'C', 'HH': 'H', 'HD13': 'H', etc.. }
    This is necessary for finding the radius of each element, when only the identifier
    of the specific position inside the chain is known.
    '''
    s_dict = {}

    with open(filepath) as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM'):
                specific = line[12:16].strip()
                _element = line[76:78].strip()
                if not (specific in s_dict):
                    s_dict[specific] = _element

    return s_dict


def vcols_from_nearest_fi(nstr, surface_obj_name, fp):

    '''
    fp = filepath for pdb
    nstr = modelremark = modelshortname = last 4 chars of pdbname

    [x] step : first store atoms as {element_name: [co,..], }
    [x] step : generate singular ordered mesh by adding vertices
               in clumps of element types. (H,H,H,H,H,H,O,O,O,O,O,O..)
    [x] step : track start and end index for each element into mapper_obj
    [x] step : get surface mesh.
    [x] step : for every vertex on surface mesh find closest
               vertex in proxy_ob, and fi value
    [ ] step : for those surface vertices which have no elements within proximity
               radius decide how to colour the mesh at that point.

    '''
    # options.
    use_nearest_after_failed_proximity = True
    normalize_and_make_numeric = True  # affects idx to fi

    ## storage
    atom_to_fi = defaultdict(list)
    proxy_obj = defaultdict(list)
    idx_to_fi = []
    mapper_obj_invert = {}
    verts = []
    surface_verts_fi = {}

    ## aliasing
    objs = bpy.data.objects
    texts = bpy.data.texts
    scene = bpy.context.scene
    meshes = bpy.data.meshes
    children = objs[nstr].children

    sd = get_specific_dict(filepath=fp)

    def generate_or_update(verts):
        # -- get or create mesh
        if "proxy_mesh" in meshes:
            mesh = meshes["proxy_mesh"]
        else:
            mesh = meshes.new("proxy_mesh")

        # -- inject mesh with verts
        bm = bmesh_from_pyverts(verts)
        bm.to_mesh(mesh)
        bm.free()

        # -- create or update object with new mesh data
        if "proxy_obj" in objs:
            obj = objs['proxy_obj']
            obj.data = mesh
        else:
            obj = objs.new("proxy_obj", mesh)
            scene.objects.link(obj)

    if not verts:
        # fill `proxy_obj` and `atom_fo_fi`
        for o in children:
            _name = o.BBInfo[12:16].strip()   # name in the chain
            _amino = o.BBInfo[17:20].strip()  # name of amino

            fi = values_fi[_amino][_name]
            co = o.location[:]

            proxy_obj[_name].append(co)
            atom_to_fi[_name].append(fi)

        # fills `idx_to_fi` and `mapper_obj`
        idx = 0
        for key in sorted(proxy_obj.keys()):
            start = len(verts)
            verts.extend(proxy_obj[key])
            end = len(verts)-1

            mapper_obj_invert[(start, end)] = key
            # mapper_obj_invert will fold range and specific element {(0, 97): 'C' ...}

            idx_to_fi.extend(atom_to_fi[key])

    if normalize_and_make_numeric:
        n = lambda fi: round(((float(fi) + 3) * 0.25), 4)
        idx_to_fi = [n(fi) for fi in idx_to_fi]

    if verts:
        generate_or_update(verts)
    else:
        print('no verts! - ending early')
        return

    if surface_obj_name and (surface_obj_name in objs):

        obj = objs.get(surface_obj_name)
        vcol_layer = get_vcol_layer(obj)

        kd = build_ktree(verts)
        max_dist = 2.7  # fudge, could test largest vdw in cloud.

        def from_closest(vidx, coordinate, mdist):
            # Each surface vertex is present in at least 3 polygons.
            # To avoid repeating proximity search with kdtree this
            # algorithm memoizes the index and resulting fi value.
            # - this assumes that a lookup in a hashtable is more
            #   efficient than kdtree + tests..for a second and third
            #   time the vertex appears in a different polygon.

            fi = surface_verts_fi.get(vidx)
            if fi:
                # return early, use a cached value
                return fi

            # -- if arrives here the value is not yet known.
            closest_elements = []

            # -- find_range returns in order of closest -> furthest
            for (co, index, dist) in kd.find_range(coordinate, mdist):
                # this is not efficient, but for testing will suffice.
                # --- translate index to element
                # k:        index ranges..(0, 20), (21, 50), etc
                # element:  the element associated that range
                for k, element in mapper_obj_invert.items():
                    if k[0] <= index <= k[1]:
                        radius = scale_vdw.get(element)
                        if not radius:
                            # element is too specific, lookup in conversion table
                            elementp = sd[element]
                            radius = scale_vdw.get(elementp)

                        radius = radius[0]
                        closest_elements.append([dist-radius, index])

                        # no need to search further
                        break

            # -- establish fi for this surface vertex
            if not closest_elements:
                # search failure..
                print('nothing within {0} of indexed surface vertex {1}'.format(mdist, vidx))
                if use_nearest_after_failed_proximity:
                    co, idx, dist = kd.find(coordinate)
                    print('found at {0}'.format(dist))
                    fi = idx_to_fi[idx]
                else:
                    fi = None
            else:
                if len(closest_elements) == 1:
                    idx = closest_elements[0][1]
                else:
                    # arrange nearest elements in order or dist-radius
                    # does not deal with interpolation of multiple close
                    list_member = sorted(closest_elements)[0]
                    idx = list_member[1]

                # at this point we have a closest idx, we know its fi.
                fi = idx_to_fi[idx]

            surface_verts_fi[vidx] = fi
            return fi

        # mesh data of surface
        surface_mesh = obj.data
        i = 0
        for poly in surface_mesh.polygons:
            for idx, vidx in zip(poly.loop_indices, poly.vertices[:]):
                coordinate = surface_mesh.vertices[vidx].co
                c = from_closest(vidx, coordinate, mdist=max_dist)
                rgb = (c, c, c) if isinstance(c, float) else (.5, .5, 1.)
                vcol_layer.data[i].color = rgb
                i += 1

        # set active to vcol fi_cols
        obj.data.vertex_colors.active = vcol_layer
