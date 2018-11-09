#!/usr/bin/env bash
echo """
nosetests
pylama
""" > .git/hooks/pre-commit
chmod +x .git/hooks/pre-commit
