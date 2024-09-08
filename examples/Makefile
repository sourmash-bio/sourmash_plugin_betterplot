.PHONY: all clean cleanall

all:
	snakemake -c 4 -p

clean:
	snakemake -c 4 --delete-all-output
	rm -f *.csv *.png

cleanall: clean all
