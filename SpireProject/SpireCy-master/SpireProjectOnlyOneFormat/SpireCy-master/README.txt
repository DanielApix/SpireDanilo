necessario modulo cython in python(2.7)
per installare: pip install Cython

per compilare
python setup.py build_ext --inplace

per avviare programma sequenziale
python run_seq.py

per avviare server programma distribuito
python run_server.py

per avviare client programma distribuito
python run_client.py

inserire percorso corretto directory che contiene le directory con i -file fasta in server.pyx e client.pyx nella directory cy_files

per macchine diverse inserire ip macchina server in file client.pyx


