.PHONY: build smoke clean

build:
	./build.sh

smoke: build
	python data/summarize.py --input results/e11_bench_summary.tsv

clean:
	rm -rf build
