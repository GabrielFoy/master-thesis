# Makefile pour compiler le programme
CXX = g++

# Flags pour le compilateur
CXXFLAGS = -Wall -Wextra -std=c++20 -I/usr/include/eigen3 -I/usr/include

# Options spécifiques à GSL
LDFLAGS = -L/usr/lib -lgsl -lgslcblas -lm -pthread

# Les fichiers sources et les fichiers objets
SRCS = main.cpp simulation.cpp visualisation_convergence.cpp
OBJS = $(SRCS:.cpp=.o)

# Le nom de l'exécutable
TARGET = main.exe visualisation_convergence.exe

# Règle par défaut
all: $(TARGET)

# Règle pour l'édition de liens pour main.exe
main.exe: main.o simulation.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# Règle pour l'édition de liens pour visualisation_convergence.exe
visualisation_convergence.exe: visualisation_convergence.o simulation.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# Règle pour compiler les fichiers objets
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Nettoyer les fichiers objets et les exécutables
clean:
	rm -f $(OBJS) $(TARGETS)