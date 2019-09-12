#version 3.7;

global_settings{assumed_gamma 1.0}
background {rgb 0}
#declare pos = <2.,1.2,1.3>;
camera {right x*1 location pos look_at <.5,.3,.5>}
light_source {pos*1.5 color rgb 1.}

#include "rotation.inc"

//Bounding Box
#declare BoundsBoxRadius = 0.002;
#declare p1=0;
#declare p2=<0,0,1>;
#declare p3=<1,0,1>;
#declare p4=<1,0,0>;
#declare p5=<0,1,0>;
#declare p6=<0,1,1>;
#declare p7=1;
#declare p8=<1,1,0>;

union{
	union{
		cylinder{p1,p2 BoundsBoxRadius}
		cylinder{p2,p3 BoundsBoxRadius}
		cylinder{p3,p4 BoundsBoxRadius}
		cylinder{p4,p1 BoundsBoxRadius}
		cylinder{p5,p6 BoundsBoxRadius}
		cylinder{p6,p7 BoundsBoxRadius}
		cylinder{p7,p8 BoundsBoxRadius}
		cylinder{p8,p5 BoundsBoxRadius}
		cylinder{p1,p5 BoundsBoxRadius}
		cylinder{p2,p6 BoundsBoxRadius}
		cylinder{p3,p7 BoundsBoxRadius}
		cylinder{p4,p8 BoundsBoxRadius}
		sphere{p1, BoundsBoxRadius}
		sphere{p2, BoundsBoxRadius}
		sphere{p3, BoundsBoxRadius}
		sphere{p4, BoundsBoxRadius}
		sphere{p5, BoundsBoxRadius}
		sphere{p6, BoundsBoxRadius}
		sphere{p7, BoundsBoxRadius}
		sphere{p8, BoundsBoxRadius}

		texture{pigment{ color rgb 1}} no_shadow
	}

	#declare R = 0.0002;
	#debug concat(concat("read in file 'data_",str(clock,1,0)),".inc'\n")


	union{
	#include concat(concat("data_",str(clock,1,0)),".inc")
	texture{pigment{color rgb <1,1,1>}} no_shadow
	}


	translate -0.5
	rotate y*clock/final_clock*ROT
	translate +0.5
}
