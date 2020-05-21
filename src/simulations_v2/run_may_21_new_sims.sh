python run_sims.py config/asymptomatic_p.yaml base asymptomatic_p_base &
python run_sims.py config/mean_time_SyID_severe.yaml base mean_time_SyID_severe_base &
python run_sims.py config/mean_time_SyID_mild.yaml base mean_time_SyID_mild_base &
python run_sims.py config/mean_time_ID.yaml base mean_time_ID_base &
python run_sims.py config/mean_time_pre_ID.yaml base mean_time_pre_ID_base &
python run_sims.py config/mean_time_E.yaml base mean_time_E_base &

python run_sims.py config/asymptomatic_p.yaml fall asymptomatic_p_fall &
python run_sims.py config/mean_time_SyID_severe.yaml fall mean_time_SyID_severe_fall &
python run_sims.py config/mean_time_SyID_mild.yaml fall mean_time_SyID_mild_fall &
python run_sims.py config/mean_time_ID.yaml fall mean_time_ID_fall &
python run_sims.py config/mean_time_pre_ID.yaml fall mean_time_pre_ID_fall &
python run_sims.py config/mean_time_E.yaml fall mean_time_E_fall &

python run_sims.py config/testing/asymptomatic_p.yaml fall asymptomatic_p_fall_testing &
python run_sims.py config/testing/mean_time_SyID_severe.yaml fall mean_time_SyID_severe_fall_testing &
python run_sims.py config/testing/mean_time_SyID_mild.yaml fall mean_time_SyID_mild_fall_testing &
python run_sims.py config/testing/mean_time_ID.yaml fall mean_time_ID_fall_testing &
python run_sims.py config/testing/mean_time_pre_ID.yaml fall mean_time_pre_ID_fall_testing &
python run_sims.py config/testing/mean_time_E.yaml fall mean_time_E_fall_testing &

