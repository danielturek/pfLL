#!/bin/bash

git pull
git add --all
echo git commit -m $*
git push

