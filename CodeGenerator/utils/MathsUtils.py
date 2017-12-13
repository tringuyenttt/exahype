##
# @file This file is part of the ExaHyPE project.
# @author ExaHyPE Group (exahype@lists.lrz.de)
#
# @section LICENSE
#
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt


import Backend




#****************************************
#****************************************
#****** Matrix/Vector operations ********
#****************************************
#****************************************

#transpose a matrix M  
def matrixTranspose(M):
    return [[M[j][i] for j in range(len(M))] for i in range(len(M[0]))]


#A dot B  
def matrixDot(A,B):
    return [[sum([A[n][k]*B[k][m] for k in range(len(B))]) for m in range(len(B[0]))] for n in range(len(A))]


#extract matrix minor without i^th row and j^th column
def matrixMinor(M,i,j):
    return [row[:j] + row[j+1:] for row in (M[:i]+M[i+1:])]    


#compute matrix determinant    
def matrixDeterminant(M):
    if len(M) == 4: # apply Leibniz formula using permutation group S_4
        return (  M[0][0]*M[1][1]*M[2][2]*M[3][3] # id
                - M[0][0]*M[1][1]*M[2][3]*M[3][2] # (23) 
                - M[0][0]*M[1][2]*M[2][1]*M[3][3] # (12) 
                + M[0][0]*M[1][2]*M[2][3]*M[3][1] # (123)
                + M[0][0]*M[1][3]*M[2][1]*M[3][2] # (132)
                - M[0][0]*M[1][3]*M[2][2]*M[3][1] # (13)
                - M[0][1]*M[1][0]*M[2][2]*M[3][3] # (01)
                + M[0][1]*M[1][0]*M[2][3]*M[3][2] # (01)(23)
                + M[0][1]*M[1][2]*M[2][0]*M[3][3] # (012)
                - M[0][1]*M[1][2]*M[2][3]*M[3][0] # (0123) = (23)o(12)o(01)
                - M[0][1]*M[1][3]*M[2][0]*M[3][2] # (0132) = (23)o(13)o(01)
                + M[0][1]*M[1][3]*M[2][2]*M[3][0] # (013)
                + M[0][2]*M[1][0]*M[2][1]*M[3][3] # (021)
                - M[0][2]*M[1][0]*M[2][3]*M[3][1] # (0231) = (13)o(23)o(02)
                - M[0][2]*M[1][1]*M[2][0]*M[3][3] # (02)
                + M[0][2]*M[1][1]*M[2][3]*M[3][0] # (023)
                + M[0][2]*M[1][3]*M[2][0]*M[3][1] # (02)(13)
                - M[0][2]*M[1][3]*M[2][1]*M[3][0] # (0213) = (13)o(12)o(02)
                - M[0][3]*M[1][0]*M[2][1]*M[3][2] # (0321) = (12)o(23)o(03)
                + M[0][3]*M[1][0]*M[2][2]*M[3][1] # (031)
                + M[0][3]*M[1][1]*M[2][0]*M[3][2] # (032)
                - M[0][3]*M[1][1]*M[2][2]*M[3][0] # (03)
                - M[0][3]*M[1][2]*M[2][0]*M[3][1] # (0312) = (12)o(13)o(03)
                + M[0][3]*M[1][2]*M[2][1]*M[3][0] # (03)(12)
               )
    if len(M) == 3:
        return M[0][0]*M[1][1]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + M[0][2]*M[1][0]*M[2][1] - M[0][1]*M[1][0]*M[2][2] - M[0][0]*M[1][2]*M[2][1] - M[0][2]*M[1][1]*M[2][0]
    if len(M) == 2:
        return M[0][0]*M[1][1] - M[0][1]*M[1][0]
    if len(M) == 1:
        return M[0][0]
    determinant = 0   
    for j in range(len(M)):
        determinant += ((-1)**j)*M[0][j]*matrixDeterminant(matrixMinor(M,0,j))
    return determinant    


#inverse matrix M    
def matrixInverse(M):
    determinant = matrixDeterminant(M)
    cofactors = []
    for i in range(len(M)):
        cofactorRow = []
        for j in range(len(M)):
            minor = matrixMinor(M,i,j)
            cofactorRow.append(((-1)**(i+j)) * matrixDeterminant(minor))
        cofactors.append(cofactorRow)
    cofactors = matrixTranspose(cofactors)
    for i in range(len(M)):
        for j in range(len(M)):
            cofactors[i][j] = cofactors[i][j]/determinant
    return cofactors

# zero-pad a vector    
def vectorPad(v,padSize):
    if padSize <= 0:
        return v
    return v + [0. for _ in range(padSize)]


# a b c   
# d e f  
# => a b c 0 d e f 0  
def matrixPadAndFlatten_RowMajor(M, padSize):
    result = []
    for i in range(len(M)):
        result += vectorPad(M[i], padSize)
    return result


# a b c   
# d e f  
# => a d 0 b e 0 c f 0    
def matrixPadAndFlatten_ColMajor(M, padSize):
    return matrixPadAndFlatten_RowMajor(matrixTranspose(M), padSize)


 
#****************************************
#****************************************
#*********** Gauss-Legendre *************
#****************************************
#****************************************  

# return Gauss-Legendre weight, point (taken from generic GaussLegendreQuadrature.cpp)
def getGaussLegendre(nDof):
    if nDof < 1:
        raise ValueError("order must be positive")
        
    if nDof > 10:
        raise ValueError("order is currently limited to 9")
        
    if nDof == 1:
        return [1.0000000000000000], [0.5000000000000000]

    if nDof == 2:
        return  [0.5000000000000000, 0.5000000000000000], \
                [0.2113248654051871, 0.7886751345948129]

    if nDof == 3:
        return  [0.2777777777777778, 0.4444444444444444, 0.2777777777777778], \
                [0.1127016653792583, 0.5000000000000000, 0.8872983346207417]

    if nDof == 4:
        return  [0.1739274225687273, 0.3260725774312732, 0.3260725774312732, 0.1739274225687273], \
                [0.0694318442029737, 0.3300094782075719, 0.6699905217924281, 0.9305681557970262]

    if nDof == 5:
        return  [0.1184634425280948, 0.239314335249683, 0.2844444444444443, 0.239314335249683, 0.1184634425280948], \
                [0.04691007703066802, 0.2307653449471584, 0.5000000000000000, 0.7692346550528415, 0.9530899229693319]

    if nDof == 6:
        return  [0.0856622461895845, 0.1803807865240695, 0.2339569672863459, 0.2339569672863459, 0.1803807865240695, 0.0856622461895845], \
                [0.03376524289842397, 0.1693953067668678, 0.3806904069584015, 0.6193095930415985, 0.8306046932331322, 0.966234757101576]

    if nDof == 7:
        return  [0.06474248308443538, 0.1398526957446382, 0.1909150252525592, 0.2089795918367344, 0.1909150252525592, 0.1398526957446382, 0.06474248308443538], \
                [0.02544604382862076, 0.1292344072003028, 0.2970774243113014, 0.5000000000000000, 0.7029225756886985, 0.8707655927996972, 0.9745539561713792]

    if nDof == 8:
        return  [0.05061426814518821, 0.1111905172266871, 0.1568533229389437, 0.1813418916891809, 0.1813418916891809, 0.1568533229389437, 0.1111905172266871, 0.05061426814518821], \
                [0.01985507175123186, 0.1016667612931866, 0.2372337950418355, 0.4082826787521751, 0.5917173212478249, 0.7627662049581645, 0.8983332387068134, 0.9801449282487682]

    if nDof == 9:
        return  [0.04063719418078751, 0.09032408034742861, 0.1303053482014677, 0.1561735385200013, 0.1651196775006297, 0.1561735385200013, 0.1303053482014677, 0.09032408034742861, 0.04063719418078751], \
                [0.01591988024618696, 0.08198444633668212, 0.1933142836497048, 0.3378732882980955, 0.5000000000000000, 0.6621267117019045, 0.8066857163502952, 0.9180155536633179, 0.984080119753813]

    if nDof == 10:
        return  [0.03333567215434358, 0.07472567457529024, 0.1095431812579912, 0.1346333596549983, 0.1477621123573766, 0.1477621123573766, 0.1346333596549983, 0.1095431812579912, 0.07472567457529024, 0.03333567215434358], \
                [0.01304673574141413, 0.06746831665550773, 0.1602952158504878, 0.2833023029353764, 0.4255628305091844, 0.5744371694908156, 0.7166976970646236, 0.8397047841495122, 0.9325316833444923, 0.9869532642585859]



#****************************************
#****************************************
#*************** ADERDG *****************
#****************************************
#****************************************

# Code taken from:    
# .. module:: aderdg
# :platform: Unix, Windows, Mac
# :synopsis: Provides routines to compute ADER-DG basis functions and operators on the unit cube.
# .. moduleauthor:: Angelika Schwarz <angelika.schwarz@tum.de>
# :synopsis: Provides routines to compute ADER-DG basis functions and operators on the unit cube.

def BaseFunc1d(xi, xin, N):
    """
    Computes the ADER-DG basis functions and their first derivative.
    
    Args:
       xi:
          The reference element point the basis functions are evaluated at.
          Here, xi refers to the greek letter that is often used as a reference element coordinate.
       xin:
          The reference element nodes corresponding to the nodal basis functions.
       N:
          Number of nodal basis functions (=order+1).
    Returns:
       phi:
          Basis function values.
       phi_xi:
          First derivatives of the basis functions.
    """
    phi    = [1.]*N 
    phi_xi = [0.]*N
    for m in range(0,N):
        for j in range(0,N):
            if j == m:
                continue 
            phi[m] = phi[m]*(xi-xin[j])/(xin[m]-xin[j])
        for i in range(0,N):
            if i == m:
                continue
            tmp = 1.;
            for j in range(0,N):
                if j == i:
                    continue
                if j == m:
                    continue
                tmp = tmp*(xi-xin[j])/(xin[m]-xin[j])
            phi_xi[m] += tmp/(xin[m]-xin[i])
    return phi, phi_xi    

def assembleStiffnessMatrix(xGPN, wGPN, N):
    """
    Computes the (reference) element stiffness matrix for an approximation of
    order N.

    Args:
       xGPN:
          Gauss-Legendre nodes (N nodes).
       wGPN:
          Gauss-Legendre weights  (N weights).
       N:
          Number of nodal basis functions (=order+1).
    Returns:
       K_xi:
          The (reference) element stiffness matrix.
    """
    # init matrix with zero
    Kxi = [[0 for _ in range(N)] for _ in range(N)]
     
    for i in range(0,N):
        phi, phi_xi = BaseFunc1d(xGPN[i], xGPN, N)
        for k in range(0,N):
            for l in range(0,N):
                Kxi[k][l] += wGPN[i]*phi_xi[k]*phi[l] 
        
    return Kxi

def assembleK1(Kxi, xGPN, N):
    """
    Computes the difference between the reference element mass operator 
    evaluated at point xi=1.0 and the element stiffness matrix.
    
    Args:
       K_xi:
          The (reference) element stiffness matrix for a approximation of 
          order N.
       xGPN:
          Gauss-Legendre nodes (N nodes).
       N:
          Order of approximation corresponding to N+1 nodal basis functions.
    Returns:
       K1:
          <unknown>
    """
    phi1, _ = BaseFunc1d(1.0, xGPN, N)
    FRm = [[0 for _ in range(N)] for _ in range(N)]
    
    for k in range(0, N):
        for l in range(0, N):
            FRm[k][l] = phi1[k]*phi1[l] 
    
    return [[FRm[i][j] - Kxi[i][j] for j in range(N)] for i in range(N)]
    
    
def assembleMassMatrix(xGPN, wGPN, N):
    """
    Computes the (reference) element mass matrix for an approximation of
    order N.

    Args:
       xGPN:
          Gauss-Legendre nodes (N nodes).
       wGPN:
          N Gauss-Legendre weights (N weights).
       N:
          Number of nodal basis functions (=order+1).
    Returns:
       M_xi:
          The (reference) element mass matrix.
    """
    # init matrix with zeros
    MM = [[0 for _ in range(N)] for _ in range(N)]
    
    for i in range(0,N):
        phi, _ = BaseFunc1d(xGPN[i], xGPN, N)
        for k in range(0,N):
            for l in range(0,N):
                MM[k][l] += wGPN[i]*phi[k]*phi[l]
      
    return MM
    
    
def assembleDiscreteDerivativeOperator(MM, Kxi):
    """
    Computes some derivative values for debugging purposes.

    Args:
       MM:
          The (reference) element mass matrix for a approximation of 
          order N.
       Kxi:
          The (reference) element stiffness matrix for a approximation of 
          order N.
       
    Returns:
       dudx:
          Derivative values for debugging purposes.
    """
    dudx = matrixDot(matrixInverse(MM),matrixTranspose(Kxi))
    return dudx    
    

def assembleFineGridProjector1d(xGPN, j, N):
    """
    Transforms the degrees of freedom located on a coarse grid edge
    nodes to degrees of freedoms located on nodes of a fine grid edge.
    The difference in levels is 1.
    
    Let us denote by P the 1d fine grid projector (=1d equidistantGridProjector). The fine grid DoF 
    are computed according to:
    
    u^{fine;j}_i =  sum_{m} P^{j}_im u^{coarse}_m
    
    Args:
       xGPN:
          Gauss-Legendre nodes (N nodes).
       j:
          Index of one the three subintervals: 0,1, or 2.
       N:
          Number of nodal basis functions (=order+1).
    Returns:
       equidistantGridProjector:
          The corresponding degrees of freedom located at nodes of an equidistant grid over (0,1).
    """
    fineGridProjector1d = [[0 for _ in range(N)] for _ in range(N)]
    
    for i in range(0, N): # Eq. basis
        phi_i, _ = BaseFunc1d((xGPN[i]+j)/3.0, xGPN, N)
        for m in range(0, N): # DG basis
            fineGridProjector1d[m][i] = phi_i[m]
    return fineGridProjector1d    
  
