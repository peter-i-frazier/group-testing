python3 run_sensitivity.py $(cat $1 | grep -v \#)  $(echo "${@:2}")
