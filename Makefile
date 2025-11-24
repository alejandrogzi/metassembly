# Copyright (c) 2020 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
# Distributed under the terms of the Apache License, Version 2.0.

# author = "Alejandro Gonzales-Irribarren"
# email = "alejandrxgzi@gmail.com"
# github = "https://github.com/alejandrogzi"
# version: 0.0.7

.PHONY: all configure

all: configure

configure:
	bash assets/configure.sh
