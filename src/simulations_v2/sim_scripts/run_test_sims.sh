python run_sims.py config/contact_tracing_recall.yaml fall contact_recall_fall &
python run_sims.py config/daily_contacts.yaml fall daily_contacts_fall &
python run_sims.py config/prevalence.yaml fall prevalence_fall &

python run_sims.py config/contact_tracing_recall.yaml fall_old_severity contact_recall_fall_old_severity &
python run_sims.py config/daily_contacts.yaml fall_old_severity daily_contacts_fall_old_severity &
python run_sims.py config/prevalence.yaml fall_old_severity prevalence_fall_old_severity &

