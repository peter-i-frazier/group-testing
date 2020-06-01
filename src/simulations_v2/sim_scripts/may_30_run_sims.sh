# june, nominal, notest
python3 run_sims.py config/contact_tracing_delay.yaml -w june -a nominal --notest &
python3 run_sims.py config/daily_contacts.yaml -w june -a nominal --notest &
python3 run_sims.py config/exposed_infection_p.yaml -w june -a nominal --notest &
python3 run_sims.py config/mild_self_reporting.yaml -w june -a nominal --notest &
python3 run_sims.py config/detailed_asymptomatic_p.yaml -w june -a nominal --notest &

# fall, nominal, notest
python3 run_sims.py config/contact_tracing_delay.yaml -w fall -a nominal --notest &
python3 run_sims.py config/daily_contacts.yaml -w fall -a nominal --notest &
python3 run_sims.py config/exposed_infection_p.yaml -w fall -a nominal --notest &
python3 run_sims.py config/mild_self_reporting.yaml -w fall -a nominal --notest &
python3 run_sims.py config/detailed_asymptomatic_p.yaml -w fall -a nominal --notest &

# fall, nominal
python3 run_sims.py config/contact_tracing_delay.yaml -w fall -a nominal &
python3 run_sims.py config/daily_contacts.yaml -w fall -a nominal &
python3 run_sims.py config/exposed_infection_p.yaml -w fall -a nominal &
python3 run_sims.py config/mild_self_reporting.yaml -w fall -a nominal &
python3 run_sims.py config/detailed_asymptomatic_p.yaml -w fall -a nominal &

# testing sims
python3 run_sims.py config/testing/testing_fraction.yaml -w fall -a nominal &
python3 run_sims.py config/testing/testing_qfnr.yaml -w fall -a nominal &
python3 run_sims.py config/testing/testing_fraction.yaml -w fall -a optimistic &
python3 run_sims.py config/testing/testing_qfnr.yaml -w fall -a optimistic &
python3 run_sims.py config/testing/testing_fraction.yaml -w fall -a pessimistic &
python3 run_sims.py config/testing/testing_qfnr.yaml -w fall -a pessimistic &

# sims added later
python3 run_sims.py config/detailed_prevalence.yaml -w june -a nominal --notest &
python3 run_sims.py config/detailed_prevalence.yaml -w fall -a nominal --notest &
python3 run_sims.py config/detailed_prevalence.yaml -w fall -a nominal &

python3 run_sims.py config/detailed_severe_self_reporting.yaml -w june -a nominal --notest &
python3 run_sims.py config/detailed_severe_self_reporting.yaml -w fall -a nominal --notest &
python3 run_sims.py config/detailed_severe_self_reporting.yaml -w fall -a nominal &

python3 run_sims.py config/detailed_contact_isolations.yaml -w june -a nominal --notest &
python3 run_sims.py config/detailed_contact_isolations.yaml -w fall -a nominal --notest &
python3 run_sims.py config/detailed_contact_isolations.yaml -w fall -a nominal &

# these are for the histograms
python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a nominal &
python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a nominal --notest &

python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a optimistic &
python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a optimistic --notest &

python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a pessimistic &
python3 run_sims.py config/single_run_via_mild_self_reporting.yaml -w fall -a pessimistic --notest &
