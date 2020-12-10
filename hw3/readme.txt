Assignment #3: Ray tracing

FULL NAME: Nikhil Johny Karuthedath
USC ID: 2900797907


MANDATORY FEATURES
------------------

<Under "Status" please indicate whether it has been implemented and is
functioning correctly.  If not, please explain the current status.>

Feature:                                 Status: finish? (yes/no)
-------------------------------------    -------------------------
1) Ray tracing triangles                  YES

2) Ray tracing sphere                     YES

3) Triangle Phong Shading                 YES

4) Sphere Phong Shading                   YES

5) Shadows rays                           YES

6) Still images                           YES (/stillImages)
   
7) Extra Credit (up to 20 points)
	Implemented antialiasing
		RAY_SAMPLES = 1, no antialiasing (Default)
		RAY_SAMPLES > 1, shoot multiple sample rays with a delta from other sample rays and then average them out 
	Implemented soft shadows
		enableSoftShadows = false, no soft shadows (Default)
		enableSoftShadows = true
		EXTRA_LIGHTS > 2, add multiple random lights with total intensity equal to 1

	Parameters used for still images:
		/stillImages/core uses RAY_SAMPLES = 1 and enableSoftShadows = false
		/stillImages/extraCredit uses RAY_SAMPLES = 2 and enableSoftShadows = true, EXTRA_LIGHTS = 4 for SIGGRAPH.scene
		/stillImages/extraCredit uses RAY_SAMPLES = 4 and enableSoftShadows = true, EXTRA_LIGHTS = 256 for every other scene