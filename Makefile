all: build/spinteny-compiled.js

build/spinteny-compiled.js: js/spinteny.js build/glow-files
	java -client -jar third-party/closure-compiler/compiler.jar --js $(shell python closure/closure/bin/build/closurebuilder.py --root=closure/ --root=js/ --namespace="Spinteny" | cat - build/glow-files) --compilation_level=ADVANCED_OPTIMIZATIONS --js_output_file $@
