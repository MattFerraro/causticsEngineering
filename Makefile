BUILD_DIR := $(PWD)/build

build:
	@docker run -it --rm -v $(PWD):/src -v $(BUILD_DIR):/root/.julia/  -w /src julia julia ./build.jl

run: build
	@docker run -it --rm -v $(PWD):/src -v $(BUILD_DIR):/root/.julia/  -w /src julia julia ./run.jl

clean:
	rm -r $(BUILD_DIR)
