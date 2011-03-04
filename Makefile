

NIM = nimrod


all: vector3D

vector3D: vector3D.nim
	${NIM} compile vector3D.nim


clean:
	rm -rf nimcache vector3D

        
.PHONY : clean
