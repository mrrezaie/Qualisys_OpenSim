# %%
from scipy.io import loadmat
import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.signal import butter, filtfilt
import os
from os.path import join
import sys
from xml.etree import ElementTree as ET
import opensim as osim
from time import sleep
from tqdm import tqdm
from shutil import rmtree
# from subprocess import call
# import matplotlib.pyplot as plt

# import winsound
# duration = 100  # milliseconds
# freq = 440  # Hz
# winsound.Beep(freq, duration)

# Level Off = 6,
# Level Critical = 5,
# Level Error = 4,
# Level Warn = 3,
# Level Info = 2,
# Level Debug = 1,
# Level Trace = 0
osim.Logger.setLevel(2)

# # add opensim dll files directory to system path
# os.environ["PATH"] += os.pathsep + 'F:\\OpenSim 4.4\\bin'


####################################################################################################
#  FUNCTIONS
####################################################################################################

def readMAT():
	# read MAT file exported by QTM
	# return marker and force dicts

	qtm = loadmat(join(WORKpath,filename), appendmat=True, simplify_cells=True)
	f = list(qtm.keys())
	qtm = qtm[f[-1]]
	f = list(qtm.keys())

	# marker data
	if 'Trajectories' in f:
		marker = dict()
		marker['labels'] = qtm['Trajectories']['Labeled']['Labels'].tolist()
		marker['used'] = int(qtm['Trajectories']['Labeled']['Count'])
		marker['start'] = 1 # int(qtm['StartFrame'])
		marker['samples'] = int(qtm['Frames'])
		marker['rate'] = int(qtm['FrameRate'])
		marker['time'] = np.linspace(0, (marker['samples']-1)/marker['rate'], marker['samples'])
		marker['frame'] = np.arange(1,marker['samples']+1)
		marker['data'] = np.transpose(qtm['Trajectories']['Labeled']['Data'][:,:3,:], (0,2,1)) # 3D data [markers, samples, x/y/z]

	# force data
	if 'Force' in f:
		force = dict()
		force['used'] = len(qtm['Force'])
		force['data'], force['location'], force['offset'] = dict(), dict(), dict()
		force['labels'] = list()
		for i in range(force['used']):
			force['samples'] = int(qtm['Force'][i]['NrOfSamples'])
			force['rate'] = int(qtm['Force'][i]['Frequency'])
			force['time'] = np.linspace(0, (force['samples']-1)/force['rate'], force['samples'])
			force['data'][i] = np.dstack(( 
										-1 * qtm['Force'][i]['Force'],
										0.001*qtm['Force'][i]['COP'], 
										-0.001*qtm['Force'][i]['Moment'])).T # 3D data [f/c/m, samples, x/y/z]
			force['location'][i] = qtm['Force'][i]['ForcePlateLocation']
			force['offset'][i] = qtm['Force'][i]['ForcePlateOffset']
			# force['labels'][i] = list()
			[force['labels'].append(f'ground_force_{i+1}_v{j}') for j in ['x', 'y', 'z']]
			[force['labels'].append(f'ground_force_{i+1}_p{j}') for j in ['x', 'y', 'z']]
			[force['labels'].append(f'ground_moment_{i+1}_m{j}') for j in ['x', 'y', 'z']]

	return marker, force




def osimCoordination(data, coordination):
    '''
    rotate 3D matrix
    data:           3D matrix e.g. [markers,row,3]
    coordination:   (e.g -xzy)
    return nothing, just rotates and update the input
    '''
    coordination = coordination.upper()
    # define which axes (X,Y,Z in order (first or second)) and how much must be rotated 
    # for the data coordinations, rot=[X,Y,Z]
    if   'XYZ' in coordination: rot=[0,0,0]; Mag1=0       ; Mag2=0 # OpenSim
    elif 'XZY' in coordination: rot=[1,0,0]; Mag1=-np.pi/2; Mag2=0
    elif 'ZYX' in coordination: rot=[0,1,0]; Mag1=+np.pi/2; Mag2=0
    elif 'ZXY' in coordination: rot=[1,0,2]; Mag1=+np.pi/2; Mag2=+np.pi/2
    elif 'YXZ' in coordination: rot=[0,2,1]; Mag1=+np.pi/2; Mag2=-np.pi
    elif 'YZX' in coordination: rot=[1,2,0]; Mag1=-np.pi/2; Mag2=-np.pi/2
    else: raise RuntimeError(f'The {coordination} is not defined properly.')
    if coordination.startswith('-'): Mag3=np.pi # for the direction of movement
    else: Mag3=0

    # function for data rotation
    if type(data)==np.ndarray:
        def rotation(Vec3, magnitude, axis):
            if   axis==0: vec = magnitude*np.array([1,0,0]) # X
            elif axis==1: vec = magnitude*np.array([0,1,0]) # Y
            elif axis==2: vec = magnitude*np.array([0,0,1]) # Z
            rotVec = R.from_rotvec(vec)
            for i in range(len(Vec3)):
                Vec3[i] = rotVec.apply(Vec3[i])

    # apply rotation function
    if Mag1!=0: rotation(data, Mag1, rot.index(1)) # Mag1
    if Mag2!=0: rotation(data, Mag2, rot.index(2)) # Mag2
    if Mag3!=0: rotation(data, Mag3, 1) # Y (vertical axis)




def rotMat(v1, v2, recompute, seq):
    '''
    calculate rotation matrix from distance vectors
    v1:             1st distance vector (2D array [row,3])
    v2:             2nd distance vector (2D array [row,3])
    recompute:      which on of the above vectors must be recomputed? 1,2,0
    seq:            sequence of the vectors (e.g. 'XYZ' for (v1,v2,v3))
    return:         rotation matrix (3D array [row,3,3])
    '''
    u1 = v1 / np.linalg.norm(v1, ord=2, axis=1).reshape(-1,1)
    u2 = v2 / np.linalg.norm(v2, ord=2, axis=1).reshape(-1,1)
    v3 = np.cross(u1, u2)
    u3 = v3 / np.linalg.norm(v3, ord=2, axis=1).reshape(-1,1)

    if recompute==1:
        u1 = np.cross(u2,u3)
        u1 = u1 / np.linalg.norm(u1, ord=2, axis=1).reshape(-1,1)
    elif recompute==2:
        u2 = np.cross(u3,u1)
        u2 = u2 / np.linalg.norm(u2, ord=2, axis=1).reshape(-1,1)
    elif recompute==0:
        pass

    seq = seq.upper()
    a = dict()
    a[0],a[1],a[2] = u1,u2,u3

    # create 3*3 orientation matrix at each time stamp
    # assign each axis to columns
    return np.dstack((a[seq.index('X')], a[seq.index('Y')], a[seq.index('Z')]))





def HarringtonHJC(RASIS, LASIS, RPSIS, LPSIS):
    '''
    calculate hipe joint centers based on regression by Harrignton
    input:      3D positions of ASIS and PSIS markers (2D array [row,3])
    
    ref:
    https://doi.org/10.1016/j.jbiomech.2006.02.003
    https://doi.org/10.1016/j.gaitpost.2015.07.004
    '''
    row = np.shape(RASIS)[0]
    midPSIS = (RPSIS + LPSIS)/2
    midASIS = (RASIS + LASIS)/2 # pelvis origin
    PW = RASIS - LASIS # z
    PD = midASIS - midPSIS # x
    pelvisTransform = rotMat(PW, PD, 2, 'ZXY') # [row,3,3] shape

    PD_norm = np.linalg.norm(PD, ord=2, axis=1)
    PW_norm = np.linalg.norm(PW, ord=2, axis=1)

    # Harrington linear regression equation
    x = -0.24 * PD_norm.reshape(-1,1) - 9.9
    y = -0.30 * PW_norm.reshape(-1,1) - 10.9
    z =  0.33 * PW_norm.reshape(-1,1) + 7.3
    # # modified Harrington using PW only by Sangeux
    # x = -0.138 * PW_norm.reshape(-1,1) - 10.4
    # y = -0.305 * PW_norm.reshape(-1,1) - 10.9
    # z =  0.33  * PW_norm.reshape(-1,1) + 7.3

    diff_R = np.hstack((x, y, z)) # [row,3] shape
    diff_L = np.hstack((x, y, -z))

    RHJC = np.empty((row,3))
    LHJC = np.empty((row,3)) 
    for i in range(row):
        RHJC[i,:] = np.dot(pelvisTransform[i,:,:], diff_R[i,:]) + midASIS[i,:]
        LHJC[i,:] = np.dot(pelvisTransform[i,:,:], diff_L[i,:]) + midASIS[i,:]
        
    return RHJC, LHJC






def find_ranges(array, thr=1):
    '''
    find the indeces (first and last) of any sequential order numbers
    thr:  a threshold for length of the indeces
    '''
    sequences = np.split(array, np.array(np.where(np.diff(array) > 1)[0]) + 1)
    l = list()
    for i in sequences:
        if len(i) > thr: l.append([np.min(i), np.max(i)])
        else: pass
    return l






def motions():
	os.system("title " + f'{filename}: Generate OpenSim motion files')

	marker, force = readMAT()

	# Process data
	if 'static' in name:

		# get specific markers
		RASIS = marker['labels'].index('RASI')
		LASIS = marker['labels'].index('LASI')
		RPSIS = marker['labels'].index('RPSI')
		LPSIS = marker['labels'].index('LPSI')
		RLFC = marker['labels'].index('RLFC')
		RMFC = marker['labels'].index('RMFC')
		LLFC = marker['labels'].index('LLFC')
		LMFC = marker['labels'].index('LMFC')
		RLMAL = marker['labels'].index('RLMAL')
		RMMAL = marker['labels'].index('RMMAL')
		LLMAL = marker['labels'].index('LLMAL')
		LMMAL = marker['labels'].index('LMMAL')

		RASIS = marker['data'][RASIS]
		LASIS = marker['data'][LASIS]
		RPSIS = marker['data'][RPSIS]
		LPSIS = marker['data'][LPSIS]
		RLFC = marker['data'][RLFC]
		RMFC = marker['data'][RMFC]
		LLFC = marker['data'][LLFC]
		LMFC = marker['data'][LMFC]
		RLMAL = marker['data'][RLMAL]
		RMMAL = marker['data'][RMMAL]
		LLMAL = marker['data'][LLMAL]
		LMMAL = marker['data'][LMMAL]
		midASIS = np.mean((RASIS, LASIS), axis=0)
		midPSIS = np.mean((RPSIS, LPSIS), axis=0)

		# find coordinate system of the laboratory based on marker data in static trial
		staticCoord = ''
		temp = np.linalg.norm((RASIS-LASIS), axis=0, ord=2)
		Z = np.argmax(temp)
		temp = np.linalg.norm((RASIS-RLMAL), axis=0, ord=2)
		Y = np.argmax(temp)
		temp = np.linalg.norm((midASIS-midPSIS), axis=0, ord=2)
		X = np.argmax(temp)

		for i in [X,Y,Z]: # [X, Y, Z]
			if   i==0: staticCoord += 'X'
			elif i==1: staticCoord += 'Y'
			elif i==2: staticCoord += 'Z'
		if np.linalg.norm(midASIS[:,X], ord=2) < np.linalg.norm(midPSIS[:,X], ord=2):
			staticCoord = '-' + staticCoord

		# rotat static
		if staticCoord != '':
			# print(f'{filename}: Lab Coordinate System:', staticCoord)
			osimCoordination(marker['data'], staticCoord)

		# Joint center estimation
		RHJC, LHJC = HarringtonHJC(RASIS, LASIS, RPSIS, LPSIS) # Harrington2007
		RKJC = np.mean((RLFC, RMFC), axis=0)
		LKJC = np.mean((LLFC, LMFC), axis=0)
		RAJC = np.mean((RLMAL, RMMAL), axis=0)
		LAJC = np.mean((LLMAL, LMMAL), axis=0)
		# append Joint center data to markers
		marker['labels'] += ['RHJC', 'LHJC', 'RKJC', 'LKJC', 'RAJC', 'LAJC']
		temp = np.dstack((RHJC.T, LHJC.T, RKJC.T, LKJC.T, RAJC.T, LAJC.T)).T
		marker['data'] = np.vstack((marker['data'], temp))
		marker['used'] += 6


	elif 'gait' in name: # gait trials

		# find coordinate system of the laboratory based on force data in dynamics trials
		dynamicCoord = ''
		for i in range(force['used']):
			if np.linalg.norm(force['data'][i][0], ord=2) > 100:
				temp = np.linalg.norm(force['data'][i][0], axis=0, ord=2)
				maxN = np.argmax(temp) # Y (vertical)
				minN = np.argmin(temp) # Z (Medial-Lateral)
				midN = np.where(np.logical_and(temp<max(temp), temp>min(temp)))[0][0] # X (Anterior-Posterior)
				for i in [midN, maxN, minN]: # [X, Y, Z]
					if   i==0: dynamicCoord += 'X'
					elif i==1: dynamicCoord += 'Y'
					elif i==2: dynamicCoord += 'Z'
				temp = force['data'][i][0][:,midN] # find Anterior-Posterior GRF component
				maxN = np.argmax(temp) # propulasion peak
				minN = np.argmin(temp) # braking peak
				if minN > maxN:
					dynamicCoord = '-' + dynamicCoord # rotate direction of gait
				break

		# rotate marker and force data
		if dynamicCoord != '':
			# print(f'{filename}: Lab Coordinate System:', dynamicCoord)
			osimCoordination(marker['data'], dynamicCoord)
			for i in range(force['used']):
				osimCoordination(force['data'][i], dynamicCoord)
				osimCoordination(force['location'][i], dynamicCoord)

		# get specific markers
		RCAL = marker['labels'].index('RCAL')
		LCAL = marker['labels'].index('LCAL')
		RTOE = marker['labels'].index('RTOE')
		LTOE = marker['labels'].index('LTOE')
		RMT5 = marker['labels'].index('RMT5')
		LMT5 = marker['labels'].index('LMT5')
		RCAL = marker['data'][RCAL]
		LCAL = marker['data'][LCAL]
		RTOE = marker['data'][RTOE]
		LTOE = marker['data'][LTOE]
		RMT5 = marker['data'][RMT5]
		LMT5 = marker['data'][LMT5]
		RF = np.mean((RCAL, RTOE, RMT5), axis=0) # right foot center marker
		LF = np.mean((LCAL, LTOE, LMT5), axis=0) # left foot center marker

		force['foot'] = dict()
		force['filt'] = dict()
		TOC[name] = dict() 
		TOTC[name] = list() 
		b, a = butter(6, 2*FC/force['rate'], 'lowpass', output='ba') # Butterworth low-pass filter

		for i in range(force['used']):
			force['filt'][i] = np.zeros_like(force['data'][i])
			if np.linalg.norm(force['data'][i][0], ord=2) > 100: # in case of fp contact

				# determine which foot contacts which force plate(s)
				loc = np.mean(force['location'][i], axis=0)
				Rtemp = np.linalg.norm((RF-loc), axis=1, ord=2) # distance between right foot center and fp center
				Ltemp = np.linalg.norm((LF-loc), axis=1, ord=2) # distance between left foot center and fp center
				thr = int(force['rate'] * 0.01) # 10 ms 
				if min(Rtemp)<min(Ltemp) and sum(Rtemp<200)>thr: # 200 mm distance to force plate center
					force['foot'][i] = 'R'
				elif min(Ltemp)<min(Rtemp) and sum(Ltemp<200)>thr:
					force['foot'][i] = 'L'

				# get time of foot contact
				V = force['data'][i][0][:,1] # vertical GRF component
				temp = np.where(V>7)[0] # 7 N threshold
				temp = find_ranges(temp, int(force['rate']*0.15)) # frames with more than 150 ms foot contact
				TOC[name][i] = (np.array(temp) -1 ) / force['rate'] # convert force frame to force time
				TOTC[name].append((np.array(temp[0]) -1 ) / force['rate'])

				# filter force data
				for j in range(3): # force/cop/moment
					for k in range(3): #x/y/z
						for m in range(len(temp)): # number of contacts
							force['filt'][i][j][temp[m][0]:temp[m][-1],k] = filtfilt(b,a, force['data'][i][j][temp[m][0]:temp[m][-1],k], padlen=100)

		TOTC[name] = np.sort(np.round(TOTC[name],4), axis=None)

	# write marker data to TRC file
	markerD = np.transpose(marker['data'], (1,0,2)) # 3D [samples, markers, x/y/z]
	markerD = markerD.reshape((marker['samples'], -1)) # 2D [samples, markers*x/y/z]
	markerD = np.hstack((marker['frame'].reshape((-1,1)), marker['time'].reshape((-1,1)), markerD))

	head0 = 'PathFileType\t4\t(X/Y/Z)\tmarkers\n' + 'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n' + f"{marker['rate']}\t{marker['rate']}\t{marker['samples']}\t{marker['used']}\tmm\t{marker['rate']}\t{marker['start']}\t{marker['samples']}\n"
	head3 = ['Frame#', 'Time']
	for i in marker['labels']:
		head3.append(i)
		head3.append('')
		head3.append('')
	head3 = '\t'.join(head3)+'\n'
	head4 = ['', '']
	for i in range(marker['used']):
		head4.append(f'X{i+1}')
		head4.append(f'Y{i+1}')
		head4.append(f'Z{i+1}')
	head4 = '\t'.join(head4)
	head = head0+head3+head4
	fmt = ['%i']
	for i in range(1 + 3*marker['used']):
		fmt.append('%.6f')
	np.savetxt(join(WORKpath,'output',f'{name}_marker.trc'), markerD, fmt=fmt, delimiter='\t', newline='\n', header=head, comments='')

	# write force data to MOT file
	if 'gait' in name:

		forceD = np.empty((force['samples'], 9*force['used']))
		for i in range(force['used']):
			temp = np.transpose(force['filt'][i], (1,0,2))
			forceD[:,9*i:9*i+9] = temp.reshape((force['samples'], -1))

		forceD = np.hstack((force['time'].reshape((-1,1)), forceD))
		head = f"forces\nversion=1\nnRows={force['samples']}\nnColumns={1+9*force['used']}\ninDegrees=yes\nendheader\n" + 'time\t' + '\t'.join(force['labels'])
		np.savetxt(join(WORKpath,'output',f'{name}_force.mot'), forceD, fmt='%.6f', delimiter='\t', newline='\n', header=head, comments='')

		# OpenSim external load file
		exl = osim.ExternalLoads()
		# exl.setName('ExtLoad')
		exl.setDataFileName(join(WORKpath,'output', f'{name}_force.mot'))
		for i in list(force['foot'].keys()):
			side = force['foot'][i]
			extForce = osim.ExternalForce()
			extForce.setName(f'{side}')
			extForce.set_applied_to_body(f'calcn_{side.lower()}')
			extForce.set_force_expressed_in_body('ground')
			extForce.set_point_expressed_in_body('ground')
			extForce.set_force_identifier(f'ground_force_{i+1}_v')
			extForce.set_point_identifier(f'ground_force_{i+1}_p')
			extForce.set_torque_identifier(f'ground_moment_{i+1}_')
			extForce.set_data_source_name(join(WORKpath,'output', f'{name}_force.mot'))
			exl.cloneAndAppend(extForce)
		exl.printToXML(join(WORKpath,'output', 'config', f'{name}_setup_external_load.xml'))






def Scale():
	os.system('title ' + f'{filename}: Opensim Scale')
	osim.Logger.removeFileSink()
	osim.Logger.addFileSink(join(WORKpath, 'output', 'config', f'{name}_log_scale.log'))

	scale = osim.ScaleTool(join(TEMPpath, 'setupScaleDefault.xml'))
	scale.setPathToSubject('')
	scale.setSubjectMass(mass)
	scale.setSubjectHeight(1000*height) # in mm
	scale.getGenericModelMaker().setModelFileName(join(TEMPpath, 'Rajagopal2015_passiveCal_hipAbdMoved.osim'))
	scale.getGenericModelMaker().setMarkerSetFileName(join(TEMPpath, 'setupMarkers.xml'))
	scale.getModelScaler().setMarkerFileName(join(WORKpath,'output', f'{name}_marker.trc'))
	time_range = osim.ArrayDouble()
	time_range.setitem(0, 0)
	time_range.setitem(1, 10)
	scale.getModelScaler().setTimeRange(time_range)
	scale.getMarkerPlacer().setMarkerFileName(join(WORKpath,'output', f'{name}_marker.trc'))
	scale.getMarkerPlacer().setTimeRange(time_range)
	scale.getMarkerPlacer().setOutputMotionFileName(join(WORKpath,'output', f'{name}_static_position.mot'))
	scale.getMarkerPlacer().setOutputModelFileName(join(WORKpath,'output', 'scaled.osim'))
	scale.printToXML(join(WORKpath,'output', 'config', f'{name}_setup_scale.xml'))
	scale.run()

	# with open(join(WORKpath, 'output', 'config', f'{name}_log_scale.log'), 'w') as f:
	# 	call(['opensim-cmd', 'run-tool', join(WORKpath, 'output', 'config', f'{name}_setup_scale.xml')], stdout=f)
	# print(f'{filename}: OpenSim Scale Done')


	# scale Maximum Isometric Force (Handsfield2014) and Maximum Contraction Velocity
	# from a function in Rajagopal sample simulation data (MATLAB)
	model = osim.Model(join(TEMPpath, 'Rajagopal2015_passiveCal_hipAbdMoved.osim'))
	scaled = osim.Model(join(WORKpath,'output', 'scaled.osim'))

	mMass = 75 # Rajagopal mass (kg)
	mHeight = 1.70 # Rajagopal height (m)
	sMass = mass
	sHeight = height
	mVolume = 47.05*mMass*mHeight + 1289.6
	sVolume = 47.05*sMass*sHeight + 1289.6

	for m,s in zip(model.getMuscles(), scaled.getMuscles()):
		mOFL = m.getOptimalFiberLength()
		sOFL = s.getOptimalFiberLength()
		scaleFactor = (sVolume/mVolume) / (sOFL/mOFL)
		sMIF = scaleFactor * m.getMaxIsometricForce()
		s.setMaxIsometricForce(sMIF)
		s.setMaxContractionVelocity(15)

	scaled.printToXML(join(WORKpath,'output', 'calibrated.osim'))






def IK():
	os.system('title ' + f'{filename}: OpenSim Inverse Kinematics')
	osim.Logger.removeFileSink()
	osim.Logger.addFileSink(join(WORKpath, 'output', 'config', f'{name}_log_IK.log'))

	ik = osim.InverseKinematicsTool(join(TEMPpath, 'setupIKDefault.xml'))
	ik.setName(f'{name}')
	ik.setResultsDir('')
	ik.set_model_file(join(WORKpath, 'output', 'calibrated.osim'))
	ik.setStartTime(0)
	ik.setEndTime(100)
	ik.set_output_motion_file(join(WORKpath, 'output', f'{name}_inverse_kinematics.mot'))
	ik.set_report_errors(True)
	ik.set_marker_file(join(WORKpath, 'output', f'{name}_marker.trc'))
	ik.printToXML(join(WORKpath, 'output', 'config', f'{name}_setup_IK.xml'))
	ik.run()

	# with open(join(WORKpath, 'output', 'config', f'{name}_log_IK.log'), 'w') as f:
	# 	call(['opensim-cmd', 'run-tool', join(WORKpath, 'output', 'config', f'{name}_setup_IK.xml')], stdout=f)
	# print(f'{filename}: OpenSim IK Done')







def ID():
	os.system('title ' + f'{filename}: OpenSim Inverse Dynamics')
	osim.Logger.removeFileSink()
	osim.Logger.addFileSink(join(WORKpath, 'output', 'config', f'{name}_log_ID.log'))

	model = osim.Model(join(WORKpath, 'output', 'calibrated.osim'))
	ind = osim.InverseDynamicsTool()
	# ind.setName('')
	ind.setModelFileName(join(WORKpath, 'output', 'calibrated.osim'))
	ind.setModel(model)
	ind.setStartTime(0)
	final = osim.Storage(join(WORKpath, 'output', f'{name}_inverse_kinematics.mot')).getLastTime()
	ind.setEndTime(final)
	ind.setCoordinatesFileName(join(WORKpath, 'output', f'{name}_inverse_kinematics.mot'))
	ind.setExternalLoadsFileName(join(WORKpath,'output', 'config', f'{name}_setup_external_load.xml'))
	ind.setLowpassCutoffFrequency(20)
	ind.setResultsDir(join(WORKpath,'output'))
	ind.setOutputGenForceFileName(join(WORKpath,'output', f'{name}_inverse_dynamics.sto'))
	fexclude = osim.ArrayStr()
	fexclude.set(0, 'Muscles')
	# fexclude.set(1, 'Actuators')
	ind.setExcludedForces(fexclude)
	# ind.set_joints_to_report_body_forces('All')
	# ind.set_output_body_forces_file(join(WORKpath, 'output', f'{name}__body_forces_at_joints.sto'))
	ind.printToXML(join(WORKpath, 'output', 'config', f'{name}_setup_ID.xml'))
	ind.run()

	# with open(join(WORKpath, 'output', 'config', f'{name}_log_ID.log'), 'w') as f:
	# 	call(['opensim-cmd', 'run-tool', join(WORKpath, 'output', 'config', f'{name}_setup_ID.xml')], stdout=f)
	# print(f'{filename}: OpenSim ID Done')







def SO():
	os.system('title ' + f'{filename}: OpenSim Static Optimization')
	osim.Logger.removeFileSink()
	osim.Logger.addFileSink(join(WORKpath, 'output', 'config', f'{name}_log_SO.log'))

	model = osim.Model(join(WORKpath, 'output', 'calibrated.osim'))
	model.initSystem()

	n = model.getNumCoordinates()
	c = model.getCoordinateSet()
	b = model.getBodySet()
	force = osim.ForceSet()
	force.setName('actuators')

	for i in range(n):
		# print(c.get(i).getName(), c.get(i).getMotionType())
		CA = osim.CoordinateActuator()
		CA.clone()
		CA.setCoordinate(c.get(i))
		if c.get(i).getMotionType() == 1: # rotational
			CA.setName(c.get(i).getName()+'_moment')
		elif c.get(i).getMotionType() == 2: # translational
			CA.setName(c.get(i).getName()+'_force')
		elif c.get(i).getMotionType() == 3: # couple
			CA.setName(c.get(i).getName()+'_force')
		# optimal force
		if c.get(i).getName().startswith('pelvis'):
			CA.setOptimalForce(200)
		elif c.get(i).getName().startswith('hip_flexion'):
			CA.setOptimalForce(2.5)
		elif c.get(i).getName().startswith('hip_adduction'):
			CA.setOptimalForce(5)
		elif c.get(i).getName().startswith('hip_rotation'):
			CA.setOptimalForce(10)
		elif c.get(i).getName().startswith('knee_angle'):
			CA.setOptimalForce(2.5)
		elif c.get(i).getName().startswith('ankle_angle'):
			CA.setOptimalForce(2.5)
		else:
			CA.setOptimalForce(1)
		force.append(CA)
	force.printToXML(os.path.join(WORKpath, 'output', 'config', 'setup_actuators.xml'))

	so = osim.StaticOptimization()
	so.setUseModelForceSet(True)
	so.setUseMusclePhysiology(True)
	so.setActivationExponent(2)
	so.setConvergenceCriterion(0.0001)
	so.setMaxIterations(100)
	model.addAnalysis(so)
	model.initSystem()

	analyze = osim.AnalyzeTool(model)
	analyze.setName(f'{name}')
	analyze.setModel(model)
	analyze.setModelFilename(join(WORKpath, 'output', 'calibrated.osim'))
	analyze.setResultsDir(join(WORKpath, 'output'))
	act = osim.ArrayStr()
	act.set(0, join(WORKpath, 'output', 'config', 'setup_actuators.xml'))
	analyze.setForceSetFiles(act)
	analyze.updateModelForces(model, join(WORKpath, 'output', 'config', 'setup_actuators.xml'))
	analyze.setCoordinatesFileName(join(WORKpath, 'output', f'{name}_inverse_kinematics.mot'))
	analyze.setExternalLoadsFileName(join(WORKpath,'output', 'config', f'{name}_setup_external_load.xml'))
	# final = osim.Storage(join(WORKpath, 'output', f'{name}_inverse_kinematics.mot')).getLastTime()
	analyze.setInitialTime(TOTC[name][0])
	analyze.setFinalTime(TOTC[name][-1])
	analyze.setLowpassCutoffFrequency(20)
	analyze.setReplaceForceSet(False)
	analyze.setSolveForEquilibrium(False)
	analyze.setLoadModelAndInput(True)
	# analyze.setStatesFileName('')
	analyze.getAnalysisSet().cloneAndAppend(so)
	analyze.printToXML(join(WORKpath, 'output', 'config', f'{name}_setup_SO.xml'))
	analyze.run()

	# with open(join(WORKpath, 'output', 'config', f'{name}_log_SO.log'), 'w') as f:
	# 	call(['opensim-cmd', 'run-tool', join(WORKpath, 'output', 'config', f'{name}_setup_SO.xml')], stdout=f)
	# print(f'{filename}: OpenSim SO Done')





def JR():
	os.system('title ' + f'{filename}: OpenSim Joint Reaction Analysis')
	osim.Logger.removeFileSink()
	osim.Logger.addFileSink(join(WORKpath, 'output', 'config', f'{name}_log_JR.log'))

	model = osim.Model(join(WORKpath, 'output', 'calibrated.osim'))
	model.initSystem()

	jr = osim.JointReaction(model)
	jr.setName('JointReaction')
	joint = osim.ArrayStr()
	joint.set(0, 'All')
	jr.setJointNames(joint)
	body = 	osim.ArrayStr()
	body.set(0, 'child')
	jr.setOnBody(body)
	frame = osim.ArrayStr()
	frame.set(0, 'child')
	jr.setInFrame(frame)
	jr.setForcesFileName(join(WORKpath, 'output', f'{name}_StaticOptimization_force.sto'))
	model.addAnalysis(jr)
	model.initSystem()

	analyze = osim.AnalyzeTool(model)
	analyze.setName(f'{name}')
	analyze.setModel(model)
	analyze.setModelFilename(join(WORKpath, 'output', 'calibrated.osim'))
	analyze.setResultsDir(join(WORKpath, 'output'))
	act = osim.ArrayStr()
	act.set(0, join(WORKpath, 'output', 'config', 'setup_actuators.xml'))
	analyze.setForceSetFiles(act)
	analyze.updateModelForces(model, join(WORKpath, 'output', 'config', 'setup_actuators.xml'))
	analyze.setCoordinatesFileName(join(WORKpath, 'output', f'{name}_inverse_kinematics.mot'))
	analyze.setExternalLoadsFileName(join(WORKpath,'output', 'config', f'{name}_setup_external_load.xml'))
	final = osim.Storage(join(WORKpath, 'output', f'{name}_inverse_kinematics.mot')).getLastTime()
	analyze.setInitialTime(TOTC[name][0])
	analyze.setFinalTime(TOTC[name][-1])
	analyze.setLowpassCutoffFrequency(20)
	analyze.setReplaceForceSet(False)
	analyze.setSolveForEquilibrium(False)
	analyze.setLoadModelAndInput(True)
	# analyze.setStatesFileName('')
	analyze.getAnalysisSet().cloneAndAppend(jr)
	analyze.printToXML(join(WORKpath, 'output', 'config', f'{name}_setup_JR.xml'))
	analyze.run()

	# with open(join(WORKpath, 'output', 'config', f'{name}_log_JR_.log'), 'w') as f:
	# 	call(['opensim-cmd', 'run-tool', join(WORKpath, 'output', 'config', f'{name}_setup_JR.xml')], stdout=f)
	# print(f'{filename}: OpenSim SO Done')



# %%

####################################################################################################
#  MAIN PROCESS
####################################################################################################

TEMPpath = os.path.dirname(sys.argv[0]) # template directory (where the script is located)
WORKpath = sys.argv[1].replace('__________', '') # working directory (where the MAT files are located)
print('Template directory:', TEMPpath)
print('Working directory:', WORKpath, '\n')

# add Geometry path to OpenSim to ignore warnings about geometry files
osim.ModelVisualizer.addDirToGeometrySearchPaths(join(TEMPpath, 'Geometry'))

# remove all outputs and creat output folder
try: rmtree(join(WORKpath, 'output'), ignore_errors=True)
except: pass
try: os.mkdir(join(WORKpath, 'output'))
except: pass
try: os.mkdir(join(WORKpath, 'output', 'config'))
except: pass

# read session.xml file
tree = ET.parse(join(WORKpath, 'session.xml'))
root = tree.getroot()
trials = [i.attrib for i in root.iter('Measurement')]

# set the static trial on top of the list 
for i,ii in enumerate(trials):
	if ii['Type'] == 'Static':
		trials.pop(i)
		trials.insert(0, ii)

# get parameters
for i in root.iter('Mass'): mass = float(i.text)
for i in root.iter('Height'): height = float(i.text)
for i in root.iter('GRF_cutoff_frequency'): FC = float(i.text)

# adjust file names
for i in trials:
	typ = i['Type']
	print(i['Filename'])
	i['Filename'] = i['Filename'][:-4]
	i['Outputname'] = i['Filename'].lower().replace(' ', '_')

TOC = dict() # times of contact (each foot)
TOTC = dict() # time of total contact (all feet)
print('\n')
# total = (len(trials)-1)*4 + 1


for i in tqdm(trials, total=len(trials), desc='OpenSim Analysis'): # Run OpenSim IK, ID, SO, and JR tools
	filename = i['Filename']
	name = i['Outputname']

	motions()

	if i['Type'] == 'Static':
		Scale()

	elif i['Type'] == 'Gait':
		IK()
		ID()
		SO()
		JR()


os.system('title ' + 'DONE')
print('\n\nFINISHED!!!')

sleep(100)
# WORKpath = 'D:\\Academic\\Codes\\Python\\MRR\\Qualisys_OpenSim\\Data\\001_test\\2022-09-01_shoe'
# TEMPpath = 'C:\\Users\\mrr\\Documents\\Qualisys_OpenSim\\Templates'
