TREE_PATH=./

run:
	python paramters ../bin/NBODY.input
	python hernquist.py

clean:
	rm -f *.pyc
