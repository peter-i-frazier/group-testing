python run_sims.py config/contact_tracing_delay.yaml base contact_delay_base &
python run_sims.py config/contact_tracing_recall.yaml base contact_recall_base &
python run_sims.py config/daily_contacts.yaml base daily_contacts_base &
python run_sims.py config/exposed_infection_p.yaml base exposed_infection_p_base &
python run_sims.py config/mild_self_reporting.yaml base mild_self_reporting_base &
python run_sims.py config/prevalence.yaml base prevalence_base &
python run_sims.py config/severe_self_reporting.yaml base severe_self_reporting_base &

python run_sims.py config/contact_tracing_delay.yaml fall contact_delay_fall &
python run_sims.py config/contact_tracing_recall.yaml fall contact_recall_fall &
python run_sims.py config/daily_contacts.yaml fall daily_contacts_fall &
python run_sims.py config/exposed_infection_p.yaml fall exposed_infection_p_fall &
python run_sims.py config/mild_self_reporting.yaml fall mild_self_reporting_fall &
python run_sims.py config/prevalence.yaml fall prevalence_fall &
python run_sims.py config/severe_self_reporting.yaml fall severe_self_reporting_fall &

python run_sims.py config/testing/contact_tracing_delay.yaml fall contact_delay_fall_testing &
python run_sims.py config/testing/contact_tracing_recall.yaml fall contact_recall_fall_testing &
python run_sims.py config/testing/daily_contacts.yaml fall daily_contacts_fall_testing &
python run_sims.py config/testing/exposed_infection_p.yaml fall exposed_infection_p_fall_testing &
python run_sims.py config/testing/mild_self_reporting.yaml fall mild_self_reporting_fall_testing &
python run_sims.py config/testing/prevalence.yaml fall prevalence_fall_testing &
python run_sims.py config/testing/severe_self_reporting.yaml fall severe_self_reporting_fall_testing &

python run_sims.py config/testing/testing_fraction.yaml fall testing_fraction_fall_testing &
python run_sims.py config/testing/testing_qfnr.yaml fall testing_qfnr_fall_testing &
