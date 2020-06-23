import Blender
from Blender import *
from Blender.Mathutils import *
#import math
from Blender import Armature as A

def getMouseLocalCoordinates(screen_x, screen_y, useMid = False):
	global myObject
#
# DESCRIPTION:
#
# actually the function name is completely misleading...
#
# the function returns a point and a direction vector in global coordinates.
# The point is the location of the virtual camera of the examined 3dwin (not a blender camera object!)
# The direction vector is the direction of a ray originating from the (vitual)screen at screencoords (screen_x, screen_y)
# and shooting towards this vitual camera position.
#
	for win3d in Window.GetScreenInfo(Window.Types.VIEW3D): # we search all 3dwins for the one containing the point (screen_x, screen_y) (could be the mousecoords for example)

		win_min_x, win_min_y, win_max_x, win_max_y = win3d.get('vertices') # calculate a few geometric extents for this window

		win_mid_x  = (win_max_x + win_min_x + 1.0) * 0.5
		win_mid_y  = (win_max_y + win_min_y + 1.0) * 0.5
		win_size_x = (win_max_x - win_min_x + 1.0) * 0.5
		win_size_y = (win_max_y - win_min_y + 1.0) * 0.5

		if useMid == True:
			screen_x = win_mid_x
			screen_y = win_mid_y
		if (win_max_x > screen_x > win_min_x) and (  win_max_y > screen_y > win_min_y):
			# if the given screencoords (screen_x, screen_y) are within the 3dwin we fount the right one...
			Window.QHandle(win3d.get('id'))
			# first we handle all pending events for this window (otherwise the matrices might come out wrong)


			# now we get a few matrices for our window...
				
			# sorry - i cannot explain here what they all do
			# - if you're not familiar with all those matrices take a look at an introduction to OpenGL...

			vmi  = Window.GetViewMatrix(); vmi.invert()  # the inverse viewing matrix
			pm   = Window.GetPerspMatrix() 			# the prespective matrix
			pmi  = Window.GetPerspMatrix(); pmi.invert()  # the inverted perspective matrix
			
			epsilon = 1e-3 # just a small value to account for floating point errors

			if (1.0 - epsilon < pmi[3][3] < 1.0 + epsilon):
				# pmi[3][3] is 1.0 if the 3dwin is in ortho-projection mode (toggled with numpad 5)

				# ortho mode: is a bit strange - actually there's no definite location of the camera ...
				# but the camera could be displaced anywhere along the viewing direction.

				d=Vector(list(Window.GetViewVector())+[0.0])

# all rays are parallel in ortho mode - so the direction vector is simply the viewing direction

				hms = Vector([(screen_x-win_mid_x) /win_size_x, (screen_y-win_mid_y) / win_size_y, 0.0, 1.0])

# these are the homogenious screencoords of the point (screen_x, screen_y) ranging from -1 to +1
				p=hms*pmi+1e6*d

# Finally we shift the position infinitely far away in
# the viewing direction to make sure the camera if outside the scene
# (this is actually a hack because this function
# is used in sculpt_mesh to initialize backface culling...)
			else:

# perspective mode: here everything is well defined - all rays converge at the camera's location
				
				dxy=[pm[3][3]*(((screen_x-win_min_x)/win_size_x)-1.0) - pm[3][0],
						 pm[3][3]*(((screen_y-win_min_y)/win_size_y)-1.0) - pm[3][1]]

				fp=Vector([pmi[0][0]*dxy[0]+pmi[1][0]*dxy[1],
					pmi[0][1]*dxy[0]+pmi[1][1]*dxy[1],
					pmi[0][2]*dxy[0]+pmi[1][2]*dxy[1]])

# fp is a global 3dpoint obtained from "unprojecting" the screenspace-point (screen_x, screen_y)
#- figuring out how to calculate this took me quite some time.
# The calculation of dxy and fp are simplified versions of my original code
#- so it's almost impossible to explain what's going on geometrically... sorry

				p=Vector([vmi[3][0],vmi[3][1],vmi[3][2],vmi[3][3]])

# the camera's location in global 3dcoords can be read directly from the inverted viewmatrix
				d=Vector(normalize_v3(sub_v3v3(p, fp))+[0.0])

# the direction vetor is simply the difference vector from the virtual camera's position
#to the unprojected (screenspace) point fp

			lp = p*myObject.getInverseMatrix()
			ld = d*myObject.getInverseMatrix() # normalize_v3
			lp.resize3D
			ld.resize3D

			return True, lp, ld
		else:
			return False,Vector(),Vector()

###########################################################################################
##sign of the funtion: Ray, Triangle, Point
def intersect_RayTriangle( R, no, cent ):
	u, v, n;             ## triangle vectors (Vector)
	dir, w0, w;          ## ray vectors (Vector)
	r, a, b;             ## params to calc ray-plane intersect (float)
 
	if no == 0:            ## triangle is degenerate
		return -1                 ## do not deal with this case
	
	#dir = R.P1 - R.P0             ## ray direction vector
	dir = R.y - R.x
	w0 = R.x - T.x
	a = -dot(n,w0)
	b = dot(n,dir)
	if fabs(b) < SMALL_NUM:     ## ray is parallel to triangle plane
		if a == 0:                ## ray lies in triangle plane
			return 2
		else:
			return 0             ## ray disjoint from plane
	
	## get intersect point of ray with triangle plane
	r = a / b
	if r < 0.0:                   ## ray goes away from triangle
		return 0                  ## => no intersect
	## for a segment, also test if (r > 1.0) => no intersect
	
	I = R.P0 + r * dir           ## intersect point of ray and plane
	
	## is I inside T?
	uu, uv, vv, wu, wv, D
	uu = dot(u,u)
	uv = dot(u,v)
	vv = dot(v,v)
	w = I - T.V0
	wu = dot(w,u)
	wv = dot(w,v)
	D = uv * uv - uu * vv
	
	## get and test parametric coords
	s, t
	s = (uv * wv - vv * wu) / D
	if s < 0.0 or s > 1.0:        ## I is outside T
		return 0
	t = (uv * wu - uu * wv) / D
	if t < 0.0 or (s + t) > 1.0:  ## I is outside T
		return 0
	
	return 1                      ## I is in T
 
def fast3dcore (screen_x, screen_y, v0, v1, v2, cu):

	# DESCRIPTION:
	# ------------
	# a do-it-all function (just for speed tests...)
	# needs tweaking: the backface culling is still not working correctly in all situations
	# ...and is only for convev meshes so far...
	# ----------------------------------------------------------------------------

	global selection_size
	global myObject, workingObject, workingMesh
	global displacement
	global mirrorMode, correspVerts
	global undoDict
	global mb_count
	global n_points, num_bones
	global points
	global start_point
	global last
	global lastadded, head_arm
	global oldx, oldy
	global armature, armature_point
	#print "Sono dentro a fast3dcore"
	mouseInView, lp,ld = getMouseLocalCoordinates(screen_x, screen_y)
	if mouseInView:
		pass
	else:
		return
	ld.resize3D()
	lp.resize3D()
	vi=Intersect( v0.co, v1.co, v2.co, ld, lp, 0)
	#print "Intersezione="
	#print vi
	#if (mb_count % 5) == 0:
		#mb.addMetaelem([0,vi.x, vi.y,vi.z, 3.0,1,3.0,1.0,1.0,1.0])				
	p = [vi.x, vi.y, vi.z]
	points.append(p)
		#print "Lunghezza array punti"
		#print len(points)
	if n_points == 0:
		cu.appendNurb([vi.x, vi.y,vi.z,1,1])
		head_arm = lastadded = last = start_point = Vector([vi.x,vi.y,vi.z])
		armature_point.append(Vector([vi.x,vi.y,vi.z]))
		oldx = screen_x
		oldy = screen_y
		n_points += 1
	else:
		if (((screen_x > (oldx+10)) or (screen_x < (oldx-10))) or ((screen_y > (oldy+10)) or (screen_y < (oldy-10)))):
			"""
			print "screen_x = %d" % screen_x
			print "screen_y = %d" % screen_y
			print "oldx = %d" % oldx
			print "oldy = %d" % oldy
			"""
			cu.appendPoint(0,[vi.x, vi.y,vi.z,1])
			last = Vector([vi.x, vi.y, vi.z])
			cu.update()
			oldx = screen_x
			oldy = screen_y
			"""
			print "Ultimo punto"
			print last
			print "Penultimo punto"
			print lastadded
			print "Vettore X"
			print vectx
			"""
			angrel = Mathutils.AngleBetweenVecs(last - lastadded,vectx)
			if angrel%10 == 0:
				print "Ci passo dopo if not angrel?"
				armature_point.append(lastadded)
				
				"""
				eb = A.Editbone()
				eb.roll = 10
				#eb.parent = arm.bones['Bone.003']
				eb.head = head_arm
				eb.tail = lastadded
				eb.options = [A.HINGE,]
				head_arm = lastadded
			#add the bone
				print 'bone num %d'%num_bones
				str = 'myNewBone%d'%num_bones
				arm.bones[str] = eb
				num_bones+=1
			#delete an old bone
			#del arm.bones['Bone.002']
			
				arm.update()  #save changes
				"""
			#angass = Mathutils.AngleBetweenVecs(last - start_point,vectx)
			#print "Angolo relativo = %d" % angrel
			print "Angolo relativo = %d" % angrel
			lastadded = Vector([vi.x, vi.y, vi.z])
			n_points += 1
	
	mb_count+=1
	selection = {}

	workingObject.makeDisplayList()
	workingMesh.update(1)
	Window.Redraw(Window.Types.VIEW3D)
	Blender.Redraw()
	
def main ():

	global selection_size
	global myMesh
	global displacement
	global myObject
	global undoMode, displaceMode, mirrorMode
	global wasinverted
	global cu_created  
	global draw
	global cu
	global num_curve
	global n_points
	global points           
	global head_arm
	#print "I'm in main"
	
	me = Mesh.Primitives.Plane(0.5)
	arm = Blender.Armature.New('Armature')
	arm.makeEditable()
	arm.drawType = Armature.OCTAHEDRON
	myObject = Object.New('Mesh')
	armobj = Object.New('Armature')
	myObject.link(me)
	sc = Scene.GetCurrent()
	sc.link(myObject)
	armobj.link(arm)
	#sc.link(armobj)
	
	"""
	eb = Armature.Editbone()
	eb.roll = 10
	eb.head = Vector(1,1,1)
	eb.tail = Vector(0,0,1)
	#eb.options = [Armature.HINGE, Armature.CONNECTED]
	arm.bones['myNewBone'] = eb
	eb.clearParent()
	eb2 = Armature.Editbone()
	eb2.roll = 10
	eb2.parent = arm.bones['myNewBone']
	eb2.head = Vector(0,0,1)
	eb2.tail = Vector(2,2,1)
	eb2.options = [Armature.HINGE, Armature.CONNECTED]
	arm.bones['myNewBone2'] = eb2
	#arm.drawType = A.STICK #set the draw type
	#arm.makeEditable() #enter editmode
	arm.update()
	#generating new editbone
	"""
	
	myObject.select(1)
	v0=me.verts[0]
	v1=me.verts[1]
	v2=me.verts[2]
	currentPerspMat = None
	lastPerspMat = Matrix(Window.GetPerspMatrix())
	doRedraw = False
	done = False
	id = None
	#try:
	#	myObject = Object.GetSelected()[0]
	id = Window.GetScreenInfo(Window.Types.VIEW3D)[0].get('id')	
	#except:
	#	print "exiting because no object was selected or no view3d window was available"
	#s	return
	#if myObject.getType() == 'Mesh':
	#	myMesh = myObject.getData()
	#	done = initializeFaceLists()
	#	scaleDispAndSelection()
	#	FindCorrespVerts()
	#else:
	#	return
	while not done:
		currentPerspMat = Window.GetPerspMatrix()
		#if viewChanged(lastPerspMat, currentPerspMat):
		#	done = initializeFaceLists()
		#	lastPerspMat = Matrix(currentPerspMat)
		#   lastPerspMat = currentPerspMat

		evt, val = Window.QRead()

		if evt in [Draw.MOUSEX, Draw.MOUSEY]:
			continue # speed-up: ignore mouse movement
		elif evt == Draw.ESCKEY: done = True
		elif (evt == Draw.LEFTMOUSE):				
			if (Window.GetKeyQualifiers() == 0):
				#prep_object()
				while Window.GetMouseButtons() == 1:	
					draw = 1
					if cu_created == 0:
						cu = Curve.New()
						#ncu = "Curve#%d" % num_curve
						ob = Object.New("Curve")
						ob.link(cu)
						sc = Blender.Scene.getCurrent()
						sc.link(ob)
						#arm = A.Armature('arm')
						#arm.drawType=A.STICK
						#arm.makeEditable()
						#armobj = Object.New('Armature')
						#armobj.link(arm)
						vv=Window.GetViewVector()
						cu_created += 1
						num_curve += 1
					mouse_x, mouse_y = Window.GetMouseCoords()
					
					fast3dcore(mouse_x, mouse_y,v0,v1,v2,cu)
					#Blender.Redraw()
				if Window.GetMouseButtons() == 0:
					print "Left Button up?"
					"""
					if draw == 1:
						head = points[0]
						tail = points[n_points-1]
						eb = Armature.Editbone()
						eb.roll = 10
						eb.head = Vector(head[0],head[1],head[2])
						eb.tail = Vector(tail[0],tail[1],tail[2])
						arm.makeEditable()
						arm.bones['myNewBone'] = eb
						eb.clearParent()
						eb2 = Armature.Editbone()
						eb2.roll = 10
						eb2.parent = arm.bones['myNewBone']
						eb2.head = Vector(tail[0],tail[1],tail[2])
						eb2.tail = Vector(2,2,1)
						eb2.options = [Armature.HINGE, Armature.CONNECTED]
						arm.bones['myNewBone2'] = eb2
						arm.update()
						stroke.append(points)
						draw = 0
						cu_created = 0
						n_points = 0
						points = []	
						#Blender.Redraw()
					"""
			elif(Window.GetKeyQualifiers() == 48) and val:
				sculpt_popup_menu()

#			else:
#				print "we should never get here"
#				print evt
#				print val
			else:
				print Window.GetKeyQualifiers()
				continue
		elif evt == Draw.REDRAW: Blender.Redraw(-1)
		elif evt == Draw.FKEY and val: print "fkey"		#Toggle displacement In/Out
		elif evt == Draw.MKEY and val: print "mkey"	#Toggle mirrorMode
		elif evt == Draw.UKEY and val:					#Toggle Undomode
			print "ukey"
		elif evt in [Draw.LEFTARROWKEY] and val: 			#shrink radius
			print "left arrow key"
		elif evt in [Draw.RIGHTARROWKEY] and val: 			#grow radius
			print "right arrow key"
		elif evt in [Draw.UPARROWKEY] and val:				#grow displacementHeight
			print "up arrow key"
		elif evt in [Draw.DOWNARROWKEY] and val:			#shrink displacementHeight
			print "down arrow key"
		else:
			#print "QAdd e QHandle?"
			sc = Blender.Scene.getCurrent()
			#print "Stroke = %d" % len(stroke)
			#for c in range(len(stroke)):
				#ob = Blender.Object.New("Mball","mb")
				#mb = Blender.Metaball.New()
				#for i in range(len(stroke[c])):
				#print "aggiungo metaballs a mb, step #%d" % i
				#print stroke[c][i]		
				#mb.addMetaelem([0, stroke[c][i][0],stroke[c][i][1],stroke[c][i][2], 2.0,1,2.0,1.0,1.0,1.0])

			ob = Blender.Object.New("Armature")
			arm = Blender.Armature.New()
			ob.link(arm)
			"""
			for m in range(len(armature)):
			"""
			for n in range(len(armature_point)):
				eb = A.Editbone()
				if n == 0:
					eb = A.Editbone()
					eb.roll = 10
						#eb.parent = arm.bones['Bone.003']
					eb.head = head_arm
					eb.tail = lastadded
						#eb.options = [A.HINGE,]
					head_arm = lastadded
						#add the bone
					arm.makeEditable()
					arm.bones['primo'] = eb
					arm.update()  #save changes
					
				#ob.link(mb)
			sc.link(ob)
			Window.QAdd(id, evt, val, 0)
			Window.QHandle(id)
			Blender.Redraw()

#def follower(mesh, camera):
	#orientamento oggetto
	#push
	#postmoltiplico
	#identita'
	#traslo -origine
	#scalo con la scala del follower
	#ruoto
	#posizione camera
	#viewup camera
	#direzione poiezione camera in vettore che poi inverto(metto - davanti)
	#poi di nuovo direzione proiezione in vettore e ne faccio il cross product con ul viewup
	#normalizzo il vettore risultante
	#faccio il cross product con il vettore direzione invertito
	#creo matrice
	#concateno
	#creo traslazione dall'origine
	#traslo
	#getmatrix(del follower)
"""
// Copy the follower's composite 4x4 matrix into the matrix provided.
void vtkFollower::GetMatrix(vtkMatrix4x4 *result)
{
  double *pos, *vup;
  double Rx[3], Ry[3], Rz[3], p1[3];
  vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
  int i;
  double distance;
  
  this->GetOrientation();
  this->Transform->Push();  
  this->Transform->PostMultiply();  
  this->Transform->Identity();

  // apply user defined matrix last if there is one 
  if (this->UserMatrix)
    {
    this->Transform->Concatenate(this->UserMatrix);
    }

  this->Transform->Translate(-this->Origin[0],
                             -this->Origin[1],
                             -this->Origin[2]);
  // scale
  this->Transform->Scale(this->Scale[0],
                         this->Scale[1],
                         this->Scale[2]);
  
  // rotate
  this->Transform->RotateY(this->Orientation[1]);
  this->Transform->RotateX(this->Orientation[0]);
  this->Transform->RotateZ(this->Orientation[2]);

  if (this->Camera)
    {
    // do the rotation
    // first rotate y 
    pos = this->Camera->GetPosition();
    vup = this->Camera->GetViewUp();

    if (this->Camera->GetParallelProjection())
      {
      this->Camera->GetDirectionOfProjection(Rz);
      Rz[0] = -Rz[0];
      Rz[1] = -Rz[1];
      Rz[2] = -Rz[2];
      }
    else
      {
      distance = sqrt(
        (pos[0] - this->Position[0])*(pos[0] - this->Position[0]) +
        (pos[1] - this->Position[1])*(pos[1] - this->Position[1]) +
        (pos[2] - this->Position[2])*(pos[2] - this->Position[2]));
      for (i = 0; i < 3; i++)
        {
        Rz[i] = (pos[i] - this->Position[i])/distance;
        }
      }
  
    // We cannot directly use the vup angle since it can be aligned with Rz:
    //vtkMath::Cross(vup,Rz,Rx);
    //vtkMath::Normalize(Rx);
    //vtkMath::Cross(Rz,Rx,Ry);       
    
    //instead use the view right angle:
    double dop[3], vur[3];
    this->Camera->GetDirectionOfProjection(dop);

    vtkMath::Cross(dop,vup,vur);
    vtkMath::Normalize(vur);

    vtkMath::Cross(Rz, vur, Ry);
    vtkMath::Normalize(Ry);
    vtkMath::Cross(Ry,Rz,Rx);

    matrix->Element[0][0] = Rx[0];
    matrix->Element[1][0] = Rx[1];
    matrix->Element[2][0] = Rx[2];
    matrix->Element[0][1] = Ry[0];
    matrix->Element[1][1] = Ry[1];
    matrix->Element[2][1] = Ry[2];
    matrix->Element[0][2] = Rz[0];
    matrix->Element[1][2] = Rz[1];
    matrix->Element[2][2] = Rz[2];

    this->Transform->Concatenate(matrix);
    }

  // translate to projection reference point PRP
  // this is the camera's position blasted through
  // the current matrix
  p1[0] = this->Origin[0] + this->Position[0];
  p1[1] = this->Origin[1] + this->Position[1];
  p1[2] = this->Origin[2] + this->Position[2];

  this->Transform->Translate(p1[0],p1[1],p1[2]);
  this->Transform->GetMatrix(result);
  
  matrix->Delete();
  this->Transform->Pop();  
} 

def test_arm():
	scn= Scene.GetCurrent()
	arm_ob= scn.getActiveObject()

	if not arm_ob or arm_ob.getType() != 'Armature':
		Draw.PupMenu('not an armature object')
		return
		for ob in scn.getChildren():
		ob.sel= 0
		arm_mat= arm_ob.matrixWorld
		arm_data= arm_ob.getData()
		bones= arm_data.bones.values()
		for bone in bones:
			bone_mat= bone.matrix['ARMATURESPACE']
			bone_mat_world= bone_mat*arm_mat
			ob_empty= Object.New('Empty', bone.name)
			scn.link(ob_empty)
			ob_empty.setMatrix(bone_mat_world)
			ob_empty.sel= 1

		arm_ob.sel= 1
		arm_ob.sel= 0

test_arm()
"""

#------------------------------------------------------------------------
mb_count = 0
num_curve = 0		
cu_created = 0
n_points = 0
n_stroke = 0
stroke = []
points = []
draw = 0
head_arm = last = lastadded = start_point = None
workingMesh = NMesh.New("temp")
workingObject = Blender.Object.New('Mesh')
workingObject.link(workingMesh)
vectx = Vector([1,0,0])
myObject = None

num_bones = 0
armature = []
armature_point = []
#------------------------------------------------------------------------

main()