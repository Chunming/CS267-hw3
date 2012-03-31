# You must use PrgEnv-pgi for Part (1) and PrgEnv-cray for Part (2).
# You must run both parts on Franklin - feel free to also experiment with Part (2) on Hopper.

CC = CC -O3
UPCC = cc -h upc -O
# If you change this, also change the mppwidth parameter in "job-knapsack-race" accordingly
NTHREADS = 4

# Removed knapsack-race from TARGETS
TARGETS=serial knapsack

all: $(TARGETS)

serial: serial.cpp
	$(CC) -o $@ $<

#knapsack-race: knapsack-race.upc
#	/global/homes/p/parkcs/franklin/bin/thrille_upcc -T$(NTHREADS) -thrille=racer $< -o $@

knapsack: knapsack.upc
	$(UPCC) -o $@ $<

clean:
	rm -f *.o $(TARGETS)
