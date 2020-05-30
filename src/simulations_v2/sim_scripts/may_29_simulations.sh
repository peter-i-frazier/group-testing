# using new run_sums.py modified 5/29
# {june, fall} x {nominal, optimistic, pessimistic}


# june, nominal 
python3 run_sims.py config/contact_tracing_delay.yaml -w june -a nominal &
python3 run_sims.py config/contact_tracing_isolations.yaml -w june -a nominal &
python3 run_sims.py config/daily_contacts.yaml -w june -a nominal &
python3 run_sims.py config/exposed_infection_p.yaml -w june -a nominal &
python3 run_sims.py config/mild_self_reporting.yaml -w june -a nominal &
python3 run_sims.py config/prevalence.yaml -w june -a nominal &
python3 run_sims.py config/severe_self_reporting.yaml -w june -a nominal &
python3 run_sims.py config/asymptomatic_p.yaml -w june -a nominal &

# june, optimistic
python3 run_sims.py config/contact_tracing_delay.yaml -w june -a optimistic &
python3 run_sims.py config/contact_tracing_isolations.yaml -w june -a optimistic &
python3 run_sims.py config/daily_contacts.yaml -w june -a optimistic &
python3 run_sims.py config/exposed_infection_p.yaml -w june -a optimistic &
python3 run_sims.py config/mild_self_reporting.yaml -w june -a optimistic &
python3 run_sims.py config/prevalence.yaml -w june -a optimistic &
python3 run_sims.py config/severe_self_reporting.yaml -w june -a optimistic &
python3 run_sims.py config/asymptomatic_p.yaml -w june -a optimistic &

# june, pessimistic
python3 run_sims.py config/contact_tracing_delay.yaml -w june -a pessimistic &
python3 run_sims.py config/contact_tracing_isolations.yaml -w june -a pessimistic &
python3 run_sims.py config/daily_contacts.yaml -w june -a pessimistic &
python3 run_sims.py config/exposed_infection_p.yaml -w june -a pessimistic &
python3 run_sims.py config/mild_self_reporting.yaml -w june -a pessimistic &
python3 run_sims.py config/prevalence.yaml -w june -a pessimistic &
python3 run_sims.py config/severe_self_reporting.yaml -w june -a pessimistic &
python3 run_sims.py config/asymptomatic_p.yaml -w june -a pessimistic &


# fall, nominal
python3 run_sims.py config/contact_tracing_delay.yaml -w fall -a nominal &
python3 run_sims.py config/contact_tracing_isolations.yaml -w fall -a nominal &
python3 run_sims.py config/daily_contacts.yaml -w fall -a nominal &
python3 run_sims.py config/exposed_infection_p.yaml -w fall -a nominal &
python3 run_sims.py config/mild_self_reporting.yaml -w fall -a nominal &
python3 run_sims.py config/prevalence.yaml -w fall -a nominal &
python3 run_sims.py config/severe_self_reporting.yaml -w fall -a nominal &
python3 run_sims.py config/asymptomatic_p.yaml -w fall -a nominal &

# fall, optimistic
python3 run_sims.py config/contact_tracing_delay.yaml -w fall -a optimistic &
python3 run_sims.py config/contact_tracing_isolations.yaml -w fall -a optimistic &
python3 run_sims.py config/daily_contacts.yaml -w fall -a optimistic &
python3 run_sims.py config/exposed_infection_p.yaml -w fall -a optimistic &
python3 run_sims.py config/mild_self_reporting.yaml -w fall -a optimistic &
python3 run_sims.py config/prevalence.yaml -w fall -a optimistic &
python3 run_sims.py config/severe_self_reporting.yaml -w fall -a optimistic &
python3 run_sims.py config/asymptomatic_p.yaml -w fall -a optimistic &

# fall, pessimistic
python3 run_sims.py config/contact_tracing_delay.yaml -w fall -a pessimistic &
python3 run_sims.py config/contact_tracing_isolations.yaml -w fall -a pessimistic &
python3 run_sims.py config/daily_contacts.yaml -w fall -a pessimistic &
python3 run_sims.py config/exposed_infection_p.yaml -w fall -a pessimistic &
python3 run_sims.py config/mild_self_reporting.yaml -w fall -a pessimistic &
python3 run_sims.py config/prevalence.yaml -w fall -a pessimistic &
python3 run_sims.py config/severe_self_reporting.yaml -w fall -a pessimistic &
python3 run_sims.py config/asymptomatic_p.yaml -w fall -a pessimistic &