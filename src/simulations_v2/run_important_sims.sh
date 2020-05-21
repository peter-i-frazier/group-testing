python run_sims.py config/contact_tracing_recall.yaml base contact_recall_base &
python run_sims.py config/daily_contacts.yaml base daily_contacts_base &
python run_sims.py config/prevalence.yaml base prevalence_base &

python run_sims.py config/contact_tracing_recall.yaml fall contact_recall_fall &
python run_sims.py config/daily_contacts.yaml fall daily_contacts_fall &
python run_sims.py config/prevalence.yaml fall prevalence_fall &

python run_sims.py config/testing/contact_tracing_recall.yaml fall contact_recall_fall_testing &
python run_sims.py config/testing/daily_contacts.yaml fall daily_contacts_fall_testing &
python run_sims.py config/testing/prevalence.yaml fall prevalence_fall_testing &

