import os


def compute_fingerprint_by_fasta(name_file):
	print("Compute fingerprint...")
	max_length_fact = 0
	factorization_file = open(name_file)
	fingerprint_file = open(name_file.replace("results", "fingerprint"), "w")
	txt_fingerprint = ""
	id_read = ''
	s = ' '

	while True:
		if len(s) == 0:
			break

		while True:
			s = factorization_file.readline().rstrip()
			if s == "":
				break
			if s[0] == '>':
				id_read = s
				break
			else:
				s = s.replace('"', '')
				s = s.replace("[", "")
				s = s.replace("]", "")
				list_fact = s.split()

			txt_fingerprint = txt_fingerprint + str(id_read) + " $ "

                        factorization = ""

			for i in range(len(list_fact)):

				factor = list_fact[i]

				if factor == "<<" or factor == ">>":
					txt_fingerprint = txt_fingerprint + factor
				else:
					txt_fingerprint = txt_fingerprint + str(len(factor))
					if max_length_fact < len(factor):
						max_length_fact = len(factor)

				if i == len(list_fact) - 1:
					txt_fingerprint = txt_fingerprint + "\n"
				else:
					txt_fingerprint = txt_fingerprint + ","
                                s += list_fact[i] + " "

			txt_fingerprint = txt_fingerprint.replace("<<", "0")
			txt_fingerprint = txt_fingerprint.replace(">>", "-1")
			fingerprint_file.write(txt_fingerprint + " $ " + s)
			txt_fingerprint = ""

	fingerprint_file.close()
	factorization_file.close()
	return max_length_fact


def kfingerprint(name_file, k):
	print("Compute k-fingerprint...")
	fingerprint_file = open(name_file)
	kfinger_file = open(name_file.replace("fingerprint", "k-finger"), 'w')

	while True:
		s = fingerprint_file.readline()
		if len(s) == 0 or s == "\n":
			break

		position = 0
		inner_fact = 0 # 0=fatt.principale, 1=fatt.secondaria
		fact_list = [] # lista lunghezze fattori effettivi
		binary_list = [] # lista fattorizzazioni applicate

		splitted = s.split()
		s = splitted[1].split(',')

		# creazione liste fact e binary
		for i in range(len(s)):
			if s[i] == '0':
				inner_fact = 1
				continue
			if s[i] == '-1':
				inner_fact = 0
				continue

			fact_list.append(int(s[i]))
			binary_list.append(inner_fact)
			
			
		if len(fact_list) < k: # se numero dei fattori piu' piccolo della finestra
			diff = k - len(fact_list)
			row = ' '.join(map(str, fact_list)) + " " + ' '.join(['*' for h in range(diff)])
			row += " $ "
			row += ' '.join(map(str, binary_list)) + " " + ' '.join(['*' for h in range(diff)])
			row += " $ "
			row += str(splitted[0]) + " " + str(position)
			
			kfinger_file.write(row + "\n")
			
		else: # scorrimento finestre
			for j in range(len(fact_list) - k + 1):
				k_finger = fact_list[j:j+k]
				b_finger = binary_list[j:j+k]

				row = ' '.join(map(str, k_finger))
				row += " $ "
				row += ' '.join(map(str, b_finger))
				row += " $ "
				row += str(splitted[0]) + " " + str(position)  # contiene id

				position += k_finger[0]

				kfinger_file.write(row + "\n")
	kfinger_file.close()

'''
dir_path_experiment = '/home/danilo/Scrivania/example'
list_runs = os.listdir(dir_path_experiment)
for run in list_runs:
	run_path = dir_path_experiment + "/" + run
	list_file = os.listdir(run_path)
	list_file = [file for file in list_file if file.startswith('results')]

	if len(list_file) == 0:  # or len(list_fasta) > 1:
		print("La directory deve contenere un file results")
		continue
'''

path = input("Fornisca percorso file dei risultati(tra apici \")\n")
compute_fingerprint_by_fasta(path)
#compute_fingerprint_by_fasta("/home/danilo/Scrivania/example/example1/results_example1.txt")
#kfingerprint('/home/danilo/Scrivania/example/example1/fingerprint_example1.txt', 4)
#kfingerprint('/mnt/c/Users/Antonio/Desktop/REGION/sampleErr/fingerprint/CFL_fingerprint_sampleErr.txt', 4)
#kfingerprint('/mnt/g/Spire/DATASET_BAM/NIST_NIST7086_H7AP8ADXX_CGTACTAG_2_NA12878/fingerprint/CFL_fingerprint_NIST_NIST7086_H7AP8ADXX_CGTACTAG_2_NA12878.txt', 4)
