import numpy as np
import scipy as sp
# import pmfgeneral as pmf
import sympy
import pypoman
from sympy.utilities.iterables import multiset_permutations
from scipy.spatial import HalfspaceIntersection
from scipy.spatial import ConvexHull, distance
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mpl_toolkits.mplot3d as a3
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
import itertools
import networkx as nx
from sympy import Plane, Point3D
from IPython.display import display, Latex, Math

#region 1 - Plotting polytope function1 #NOTE: Sometimes, it is not working properly. These functions are taken from stackexchange
class Faces():
    def __init__(self,tri, sig_dig=12, method="convexhull"):
        self.method=method
        self.tri = np.around(np.array(tri), sig_dig)
        self.grpinx = list(range(len(tri)))
        norms = np.around([self.norm(s) for s in self.tri], sig_dig)
        _, self.inv = np.unique(norms,return_inverse=True, axis=0)

    def norm(self,sq):
        cr = np.cross(sq[2]-sq[0],sq[1]-sq[0])
        return np.abs(cr/np.linalg.norm(cr))

    def isneighbor(self, tr1,tr2):
        a = np.concatenate((tr1,tr2), axis=0)
        return len(a) == len(np.unique(a, axis=0))+2

    def order(self, v):
        if len(v) <= 3:
            return v
        v = np.unique(v, axis=0)
        n = self.norm(v[:3])
        y = np.cross(n,v[1]-v[0])
        y = y/np.linalg.norm(y)
        c = np.dot(v, np.c_[v[1]-v[0],y])
        if self.method == "convexhull":
            h = ConvexHull(c)
            return v[h.vertices]
        else:
            mean = np.mean(c,axis=0)
            d = c-mean
            s = np.arctan2(d[:,0], d[:,1])
            return v[np.argsort(s)]

    def simplify(self):
        for i, tri1 in enumerate(self.tri):
            for j,tri2 in enumerate(self.tri):
                if j > i: 
                    if self.isneighbor(tri1,tri2) and \
                       self.inv[i]==self.inv[j]:
                        self.grpinx[j] = self.grpinx[i]
        groups = []
        for i in np.unique(self.grpinx):
            u = self.tri[self.grpinx == i]
            u = np.concatenate([d for d in u])
            u = self.order(u)
            groups.append(u)
        return groups

def plot_poly(vertices,azimuth=10,elev=10):
    hull = ConvexHull(vertices.T)
    simplices = hull.simplices

    org_triangles = [vertices.T[s] for s in simplices]
    
    f = Faces(org_triangles)
    g = f.simplify()

    ax = a3.Axes3D(plt.figure())

    colors = list(map("C{}".format, range(len(g))))

    pc = a3.art3d.Poly3DCollection(g,  facecolor=colors, 
                                       edgecolor="k", alpha=0.9)
    ax.add_collection3d(pc)

    ax.dist=10
    ax.azim=azimuth
    ax.elev=elev
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim3d([-1,1])
    plt.show()
#endregion 

#region 2 - Plotting polytope functions 2. These functions are taken from stackexchange
def simplify(triangles):
    """
    Simplify an iterable of triangles such that adjacent and coplanar triangles form a single face.
    Each triangle is a set of 3 points in 3D space.
    """

    # create a graph in which nodes represent triangles;
    # nodes are connected if the corresponding triangles are adjacent and coplanar
    G = nx.Graph()
    G.add_nodes_from(range(len(triangles)))
    for ii, a in enumerate(triangles):
        for jj, b in enumerate(triangles):
            if (ii < jj): # test relationships only in one way as adjacency and co-planarity are bijective
                if is_adjacent(a, b):
                    if is_coplanar(a, b, np.pi / 180.):
                        G.add_edge(ii,jj)

    # triangles that belong to a connected component can be combined
    components = list(nx.connected_components(G))
    simplified = [set(flatten(triangles[index] for index in component)) for component in components]

    # need to reorder nodes so that patches are plotted correctly
    reordered = [reorder(face) for face in simplified]

    return reordered

def is_adjacent(a, b):
    return len(set(a) & set(b)) == 2 # i.e. triangles share 2 points and hence a side

def is_coplanar(a, b, tolerance_in_radians=0):
    a1, a2, a3 = a
    b1, b2, b3 = b
    plane_a = Plane(Point3D(a1), Point3D(a2), Point3D(a3))
    plane_b = Plane(Point3D(b1), Point3D(b2), Point3D(b3))
    if not tolerance_in_radians: # only accept exact results
        return plane_a.is_coplanar(plane_b)
    else:
        angle = plane_a.angle_between(plane_b).evalf()
        angle %= np.pi # make sure that angle is between 0 and np.pi
        return (angle - tolerance_in_radians <= 0.) or \
            ((np.pi - angle) - tolerance_in_radians <= 0.)

flatten = lambda l: [item for sublist in l for item in sublist]

def reorder(vertices):
    """
    Reorder nodes such that the resulting path corresponds to the "hull" of the set of points.

    Note:
    -----
    Not tested on edge cases, and likely to break.
    Probably only works for convex shapes.

    """
    if len(vertices) <= 3: # just a triangle
        return vertices
    else:
        # take random vertex (here simply the first)
        reordered = [vertices.pop()]
        # get next closest vertex that is not yet reordered
        # repeat until only one vertex remains in original list
        vertices = list(vertices)
        while len(vertices) > 1:
            idx = np.argmin(get_distance(reordered[-1], vertices))
            v = vertices.pop(idx)
            reordered.append(v)
        # add remaining vertex to output
        reordered += vertices
        return reordered

def get_distance(v1, v2):
    v2 = np.array(list(v2))
    difference = v2 - v1
    ssd = np.sum(difference**2, axis=1)
    return np.sqrt(ssd)

def Plot3DPolytope(verts,dist=10,azim=50,elev=10,xmin=-1,xmax=1,ymin=-1,ymax=1,zmin=-1,zmax=1):
    hull = ConvexHull(verts)
    faces = hull.simplices

    ax = a3.Axes3D(plt.figure())
    ax.dist = dist
    ax.azim = azim
    ax.elev = elev

    ax.set_xlim3d([xmin, xmax])
    ax.set_ylim3d([ymin, ymax])
    ax.set_zlim3d([zmin, zmax])

    triangles = []
    for s in faces:
        sq = [
            (verts[s[0], 0], verts[s[0], 1], verts[s[0], 2]),
            (verts[s[1], 0], verts[s[1], 1], verts[s[1], 2]),
            (verts[s[2], 0], verts[s[2], 1], verts[s[2], 2])
        ]
        triangles.append(sq)

    new_faces = simplify(triangles)
    for sq in new_faces:
        f = a3.art3d.Poly3DCollection(list(sq),alpha=0.4)
        f.set_color('b')#colors.rgb2hex(sp.rand(3)))
        f.set_edgecolor('k')
        f.set_alpha(0.4)
        ax.add_collection3d(f)

def Plot2DPolytope(verts):
    pypoman.plot_polygon(verts)

def simpleplot3dpoly(cube,elev=10,azim = 0, dist = 10,color='b', figsize = (8,5)):
    fig = plt.figure(figsize = figsize)
    ax = fig.add_subplot(111, projection="3d")
    hull = ConvexHull(cube)
    # draw the polygons of the convex hull
    for s in hull.simplices:
        tri = Poly3DCollection(cube[s])
        tri.set_color(color)
        tri.set_alpha(0.5)
        ax.add_collection3d(tri)
    # draw the vertices
    ax.scatter(cube[:, 0], cube[:, 1], cube[:, 2], marker='o', color='purple')

    ax.view_init(elev=elev, azim=azim)
    ax.dist = dist
    plt.show()

#endregion

#region 3 - Set Theoretic Functions

def findsubsets(S,m):
    return set(itertools.combinations(S, m))

def swapped(l, p1, p2):
    r = l[:]
    r[p1], r[p2] = r[p2], r[p1]
    return r

def transpositions(A):
    return [swapped(A,i,j) for i,j in list(findsubsets(A,2))]

def elementary_perm_matrix(size,frm,to):
    P = np.identity(size)
    P[:,[frm,to]] = P[:,[to,frm]]
    return P

def check_sign_perm(perm_mat):
    n = perm_mat.shape[0]
    if (np.linalg.norm((np.abs(perm_mat)>1e-6).sum(axis = 0) - np.ones((1,n))) < 1e-6):
        return True
    else: 
        return False
#endregion

#region 4 - Identifiability Functions 
def CheckPolytopeNew(V):
    n = V.shape[0]
    m = V.shape[1]
    checkval = 0
    indxs = np.arange(m)
    print('{:d} Vertices:'.format(m))
    Vp = np.linalg.pinv(V)
    for i,indxp in enumerate(multiset_permutations(indxs)):
        Vindxp = V[:, indxp]
        G = Vindxp@Vp
        if (np.linalg.norm(G@V-Vindxp) < 1e-6):
            GG = G*(abs(G) > 1e-3)
            # if (np.linalg.norm(np.abs(np.sign(GG)).sum(axis = 0) - np.ones((1,n))) > 1e-6):
            #     # print('Not identifiable, the index is {}'.format(i))
            #     checkval += 1
            if not check_sign_perm(GG):
                checkval +=1

        
    if (checkval > 0):
        print('Not Identifiable!')
    else:
        print('Identifiable!')
        
def Check_Identifiability_Poly1(V,print_result = False):
    n = V.shape[0]
    m = V.shape[1]
    checkval = 0
    no_conclusion = 0
    indxs = np.arange(m)
    print('{:d} Vertices:'.format(m))

    Vp = np.linalg.pinv(V)
    for i,indxp in enumerate(transpositions(indxs)):
        Vindxp = V[:, indxp]
        G = Vindxp@Vp
        if (np.linalg.norm(G@V-Vindxp) < 1e-6):
            GG = G*(abs(G) > 1e-3)
            if (np.linalg.norm(np.abs(GG).sum(axis = 0) - np.ones((1,n))) > 1e-6):
                if print_result: 
                    print('Not identifiable, transposition index is {}'.format(i))
                checkval += 1
        else:
            no_conclusion += 1
            
    if (checkval > 0):
        print('Not Identifiable!')
        return False
    elif ((checkval == 0) & (no_conclusion>0)):
        print('We cannot conclude anything from this function')
        return None
    else:
        print('Identifiable!')
        return True

def Identifiability_Spoiler_Matrix(V):
    n = V.shape[0]
    m = V.shape[1]
    checkval = 0
    indxs = np.arange(m)
    print('{:d} Vertices:'.format(m))
    Vp = np.linalg.pinv(V)
    for i,indxp in enumerate(multiset_permutations(indxs)):
        Vindxp = V[:, indxp]
        G = Vindxp@Vp
        if (np.linalg.norm(G@V-Vindxp) < 1e-6):
            GG = G*(abs(G) > 1e-3)
            if (np.linalg.norm(np.abs(GG).sum(axis = 0) - np.ones((1,n))) > 1e-6):
                print('Idenfiability spoiler matrix is found')
                checkval += 1
                break
    if checkval>0:
        return GG
    else:
        print('No identifiability spoiler matrix is found')
        return None

#endregion


#region last - Rest of the functions
def display_matrix(array):
    data = ''
    for line in array:
        if len(line) == 1:
            data += ' %.3f &' % line + r' \\\n'
            continue
        for element in line:
            data += ' %.3f &' % element
        data += r' \\' + '\n'
    display(Math('\\begin{bmatrix} \n%s\end{bmatrix}' % data))

#endregion

