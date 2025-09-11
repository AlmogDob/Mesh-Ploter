CCHECKS = -fsanitize=address
CWARNINGS = -Wall -Wextra -Wuninitialized 
CFLAGS = $(CWARNINGS) -lm -lSDL2 -lSDL2_ttf -O3
O_FILES_MAIN = ./build/main.o ./build/mesher.o ./build/solver.o
O_FILES_animate_solver = ./build/animate_solver.o ./build/mesher.o ./build/solver.o

# IN_FILE=input.txt OUT_DIR=./results make main
main: build_and_link_main
	@echo
	./build/main $(IN_FILE) $(OUT_DIR)

	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_MAIN) ./build/main

	@echo
	@echo [INFO] done

build_main: 
	@echo [INFO] building main
	@gcc -c ./src/main.c $(CFLAGS) -o ./build/main.o

build_mesher: ./src/mesher.c
	@echo [INFO] building mesher
	@gcc -c ./src/mesher.c $(CFLAGS) -o ./build/mesher.o

build_solver: ./src/solver.c
	@echo [INFO] building solver
	@gcc -c ./src/solver.c $(CFLAGS) -o ./build/solver.o

link_main: $(O_FILES_MAIN)
	@echo [INFO] linking
	@gcc $(O_FILES_MAIN) $(CFLAGS) -o ./build/main

build_and_link_main: build_mesher build_solver build_main link_main

clean_main:
	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_MAIN) ./build/main


debug_main: debug_build_mesher debug_build_main debug_build_solver link_main
	gdb ./build/main

	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_MAIN)

	@echo
	@echo [INFO] done

debug_build_main: 
	@echo [INFO] building main
	@gcc -c ./src/main.c $(CFLAGS) -ggdb -o ./build/main.o

debug_build_mesher: ./src/mesher.c
	@echo [INFO] building mesher
	@gcc -c ./src/mesher.c $(CFLAGS) -ggdb -o ./build/mesher.o

debug_build_solver: ./src/solver.c
	@echo [INFO] building solver
	@gcc -c ./src/solver.c $(CFLAGS) -ggdb -o ./build/solver.o

profile_main: profile_build_mesher profile_build_solver profile_build_main profile_link_main
	./build/main $(IN_FILE) $(OUT_DIR)
	@echo
	gprof ./build/main gmon.out | /home/almog/.local/bin/gprof2dot | dot -Tpng -Gdpi=200 -o output.png
	# @sleep 0.1
	imview ./output.png
	@echo
	rm ./gmon.out ./output.png 
	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_MAIN)

profile_link_main: $(O_FILES_MAIN)
	@echo [INFO] linking
	@gcc $(O_FILES_MAIN) $(CFLAGS) -p -ggdb -o ./build/main

profile_build_main: 
	@echo [INFO] building main
	@gcc -c ./src/main.c $(CFLAGS) -p -ggdb -o ./build/main.o

profile_build_mesher: ./src/mesher.c
	@echo [INFO] building mesher
	@gcc -c ./src/mesher.c $(CFLAGS) -p -ggdb -o ./build/mesher.o

profile_build_solver: ./src/solver.c
	@echo [INFO] building solver
	@gcc -c ./src/solver.c $(CFLAGS) -p -ggdb -o ./build/solver.o

# valgrind -s --leak-check=full ./main
# cloc --exclude-lang=JSON,make .

############################################################

animate_solver: build_and_link_animate_solver
	@echo
	./build/animate_solver $(IN_FILE) $(OUT_DIR)

	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_animate_solver) ./build/animate_solver

	@echo
	@echo [INFO] done

build_animate_solver: 
	@echo [INFO] building animate_solver
	@gcc -c ./src/animate_solver.c $(CFLAGS) -o ./build/animate_solver.o

link_animate_solver: $(O_FILES_animate_solver)
	@echo [INFO] linking
	@gcc $(O_FILES_animate_solver) $(CFLAGS) -o ./build/animate_solver

build_and_link_animate_solver: build_mesher build_solver build_animate_solver link_animate_solver

clean_animate_solver:
	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_animate_solver) ./build/animate_solver


debug_animate_solver: debug_build_mesher debug_build_animate_solver debug_build_solver link_animate_solver
	gdb ./build/animate_solver

	@echo
	@echo [INFO] removing build files
	rm -r $(O_FILES_animate_solver)

	@echo
	@echo [INFO] done

debug_build_animate_solver: 
	@echo [INFO] building animate_solver
	@gcc -c ./src/animate_solver.c $(CFLAGS) -ggdb -o ./build/animate_solver.o