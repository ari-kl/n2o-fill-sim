@echo off
rem Opens n2o-fill-sim on Windows. Arguments are passed on: a scene file, --help, and so on.
julia --project="%~dp0." "%~dp0bin\n2o-fill-sim.jl" %*
