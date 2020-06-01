
python3 run_sims.py config/detailed_prevalence.yaml -w june -a nominal --notest &
python3 run_sims.py config/detailed_prevalence.yaml -w june -a nominal --notest &
python3 run_sims.py config/detailed_prevalence.yaml -w fall -a nominal &

python3 run_sims.py config/detailed_severe_self_reporting.yaml -w june -a nominal --notest &
python3 run_sims.py config/detailed_severe_self_reporting.yaml -w fall -a nominal --notest &
python3 run_sims.py config/detailed_severe_self_reporting.yaml -w fall -a nominal &

python3 run_sims.py config/detailed_contact_isolations.yaml -w june -a nominal --notest &
python3 run_sims.py config/detailed_contact_isolations.yaml -w fall -a nominal --notest &
python3 run_sims.py config/detailed_contact_isolations.yaml -w fall -a nominal &

python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a nominal &
python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a nominal --notest &

python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a optimistic &
python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a optimistic --notest &

python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a pessimistic &
python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a pessimistic --notest &
