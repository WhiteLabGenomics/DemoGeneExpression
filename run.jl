#!/usr/bin/env -S julia --project
using Pkg
Pkg.instantiate()
using DemoGeneExpression: launch
itm = launch()
sleep(1.)
run(`xdg-open "$(itm.url)"`)
wait(itm.task)
