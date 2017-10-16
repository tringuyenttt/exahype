pdi = self.GetInput()

pdo = self.GetOutput()
newPoints = vtk.vtkPoints()
   
numPoints=pdi.GetNumberOfPoints()

vals= pdi.GetPointData().GetVectors("Q")

points_dict={}
  
for point_id in range(0, numPoints):
    
    q  = vals.GetTuple(point_id) 
    x_q = q[-3]
    y_q = q[-2]
    z_q = q[-1]    
    
    newPoints.InsertPoint(point_id, x_q, y_q, z_q)
  
pdo.SetPoints(newPoints)
