UNAME_S := $(shell uname -s)


DEBUG_CFLAGS:=-g3 -O0 -pedantic -Wall -Wextra -fmessage-length=0 -std=gnu++0x `pkg-config sdl  --cflags`
RELEASE_CFLAGS:= -g -O3 -pedantic -Wall -Wextra -fmessage-length=0 -std=gnu++0x `pkg-config sdl  --cflags` -fopenmp

CLANG_CFLAGS:= -g -O3 -pedantic -Wall -Wextra -fmessage-length=0 -std=c++0x `pkg-config sdl  --cflags` -DNO_OMP 

# ifeq ($(UNAME_S),Linux)
# 	RELEASE_CFLAGS += -march=native
# endif
ifeq ($(UNAME_S),Darwin)
	RELEASE_CFLAGS += -framework OpenGL  `pkg-config --cflags opencv`
	DEBUG_CFLAGS += -framework OpenGL `pkg-config --cflags opencv`
endif

LDFLAGS:= -IInclude -LLib -lSDL_image -lSDL_mixer -lGL `pkg-config sdl  --libs`
ifeq ($(UNAME_S),Linux)
	LDFLAGS += -lv4l2
endif
ifeq ($(UNAME_S),Darwin)
	LDFLAGS += `pkg-config --libs opencv`
endif


SOURCES := src/main.cpp
SOURCES += src/CVTKTools.cpp
ifeq ($(UNAME_S),Linux)
	SOURCES += src/v4l/imageReader.cpp
endif

all:	release

install_deb_packages:
	sudo apt-get install libsdl-image1.2-dev libsdl1.2-dev libv4l-dev

kinect:
	mkdir -p ./build
	$(CXX) $(RELEASE_CFLAGS) $(SOURCES) -o build/fa_2013_release $(LDFLAGS) -DKINECT_INCLUDED -lOpenNI

parallel:
	mkdir -p ./build
	$(CXX) $(RELEASE_CFLAGS) $(SOURCES) -DPARALLEL -o build/fa_2013_release $(LDFLAGS) -ltbb

parallel_kinect:
	mkdir -p ./build
	$(CXX) $(RELEASE_CFLAGS) $(SOURCES) -DPARALLEL -o build/fa_2013_release $(LDFLAGS) -ltbb -DKINECT_INCLUDED -lOpenNI

parallel_debug:
	mkdir -p ./build
	$(CXX) $(DEBUG_CFLAGS) $(SOURCES) -DPARALLEL=1 -o build/fa_2013_debug $(LDFLAGS)

release:
	mkdir -p ./build
	$(CXX) $(RELEASE_CFLAGS) $(SOURCES) -o build/fa_2013_release $(LDFLAGS)

game2:
	mkdir -p ./build
	$(CXX) $(RELEASE_CFLAGS) $(SOURCES) -DGAME2 -o build/fa_2013_game2 $(LDFLAGS)

game2_kinect:
	mkdir -p ./build
	$(CXX) $(RELEASE_CFLAGS) $(SOURCES) -DGAME2 -o build/fa_2013_game2 -DKINECT_INCLUDED $(LDFLAGS) -lOpenNI

game2_parallel_kinect:
	mkdir -p ./build
	$(CXX) $(RELEASE_CFLAGS) $(SOURCES) -DPARALLEL -DGAME2 build/fa_2013_game2 $(LDFLAGS) -ltbb -DKINECT_INCLUDED -lOpenNI

navier:
	echo "You have to built trilinos first: run do_configure.sh and built/make"
	mkdir -p ./build 	
	$(CXX) $(RELEASE_CFLAGS) -I trilinos/built/packages/epetra/src -I trilinos/packages/triutils/src -I trilinos/packages/triutils/src -I trilinos/packages/teuchos/core/src -I trilinos/packages/teuchos/core/src -I trilinos/packages/epetra/src -I trilinos/built/packages/aztecoo/src -I trilinos/packages/aztecoo/src $(SOURCES) -DNAVIER_STOKES -o build/fa_2013_navier $(LDFLAGS) -L trilinos/built/packages/epetra/src -L trilinos/built/packages/aztecoo/src -L trilinos/built/packages/triutils/src -L trilinos/built/packages/teuchos/core/src -laztecoo -lepetra -lblas -lteuchoscore -llapack -ltriutils

navier-debug:
	echo "You have to built trilinos first: run do_configure.sh and built/make"
	mkdir -p ./build 	
	$(CXX) $(DEBUG_CFLAGS) -I trilinos/built/packages/epetra/src -I trilinos/packages/triutils/src -I trilinos/packages/triutils/src -I trilinos/packages/teuchos/core/src -I trilinos/packages/teuchos/core/src -I trilinos/packages/epetra/src -I trilinos/built/packages/aztecoo/src -I trilinos/packages/aztecoo/src $(SOURCES) -DNAVIER_STOKES -o build/fa_2013_navier $(LDFLAGS) -L trilinos/built/packages/epetra/src -L trilinos/built/packages/aztecoo/src -L trilinos/built/packages/triutils/src -L trilinos/built/packages/teuchos/core/src -laztecoo -lepetra -lblas -lteuchoscore -llapack -ltriutils


debug:
	mkdir -p ./build
	$(CXX) $(DEBUG_CFLAGS) $(SOURCES) -o build/fa_2013_debug $(LDFLAGS) 

tools:
	mkdir -p ./build
	$(CXX) $(RELEASE_CFLAGS) src/getcolor.c -o build/getcolor -lX11 

fft:
	mkdir -p ./build
	$(CXX) $(RELEASE_CFLAGS) -DUSE_FFT -lfftw3 $(SOURCES) -o build/fa_2013_release $(LDFLAGS)

clang:
	mkdir -p ./build
	clang++ $(CLANG_CFLAGS) $(SOURCES) -o build/fa_2013_release $(LDFLAGS)

clean:
	rm -rf build Debug Release
