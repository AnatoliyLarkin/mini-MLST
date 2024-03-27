# -*- coding: utf-8 -*-

#!pip install pandas
#!pip install regex
#!pip install biopython

import regex as re
import pandas as pd
import Bio
from Bio import AlignIO, Align

df = pd.read_csv('MLST_db_v2.csv', sep = ';')


def validation (request,regime):
  if request == '':
    return(False)
  elif regime == '1' and len(re.findall(r'[ATGCN]', request)) == len(request):
    return('st_from_profile')
  elif regime == '2' and request.isdigit() == True:
    return('profile_from_st')
  elif regime == '3' and request.isdigit() == True:
    return('genogroup_from_st')
  else:
    return(False)


def Profile_to_ST (request,df):
  aligner = Align.PairwiseAligner()
  aligner.match_score = 1.0
  aligner.mismatch_score = -1
  aligner.gap_score = -100

  flag = False
  for index, row in df.iterrows():
    if row['MLST_Profile'] == request:
      flag = True
      print(f'------------------------------------------------ \n')
      print(f'Requested profile : {request}')
      print('Exact match found')
      print(f"Predicted MLST: {row['MLST']}")
      alignments = aligner.align(row['MLST_Profile'], request)
      for i in alignments:
        print(i)
      print ('Match 18/18')
      print(f'------------------------------------------------ \n')

      break
  if flag == False:
    alignments = {}
    scorearray = {}
    genogroupchecker = str()
    print(f'Requested profile : {request}')
    print('Exact matches not found')
    for index, row in df.iterrows():
      if str(row['Genogroup']) != genogroupchecker or genogroupchecker == 'singleST':
        ali = aligner.align(row['MLST_Profile'], request)
        alignments[row['MLST']] = ali
        scorearray [row ['MLST']] = ali.score
        genogroupchecker = row['Genogroup']

    sorted_scores = dict(sorted(scorearray.items(), key=lambda x:x[1])[::-1])

    output = dict()
    for i in sorted_scores:
      if sorted_scores[i] >= 14:
        print(f'Similar MLST: {i}')
        print(alignments[i][0])
        print(f'Match {int((18 + sorted_scores[i])/2)}/18')
        print(f'------------------------------------------------ \n')













def ST_to_Profile(request,df):
  flag = False
  snpinfo = ''
  for index,row in df.iterrows():
    if str(row['MLST']) == str(request):
      st_data = row[:21].to_dict()
      print('Requested MLST found.')
      print('Results:')
      print(f'------------------------------------------------ \n')
      print(f"Selected MLST: {st_data['MLST']}")


      for i in st_data:
        if i != 'MLST' and i != 'MLST_Profile' and i != 'Genogroup':
          #print(f"Position {i}: SNP - {st_data[i]}")
          snpinfo += st_data[i]
      print(f"SNP profile for MLST {st_data['MLST']}" )
      print('-' * len('| 3 | 2 | 5 | 6 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 2 | 2 | 3 |'))
      print('| 3 | 2 | 5 | 6 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 2 | 2 | 3 |')
      print('| 4 | 0 | 6 | 1 | 1 | 2 | 5 | 6 | 6 | 6 | 6 | 8 | 9 | 0 | 2 | 5 | 2 | 2 |')
      print('|   | 5 | 5 | 1 | 1 | 7 | 2 | 4 | 4 | 8 | 9 | 8 | 3 | 2 | 8 | 1 | 7 | 7 |')
      print('|   |   |   |   | 2 | 1 | 6 | 6 | 9 | 5 | 1 | 6 | 7 | 8 | 2 | 3 | 5 | 6 |')
      print('-' * len('| 3 | 2 | 5 | 6 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 2 | 2 | 3 |'))

      snpinfo = list(snpinfo)
      snpinfo = '|'.join(snpinfo)
      snpinfo = f'|{snpinfo}|'
      print(*snpinfo)

      print('\n')
      print(f"MLST profile: {st_data['MLST_Profile']}")
      print('\n')
      print(f'------------------------------------------------ \n')
      flag = True
      break
  if flag == False:
    print('Requested MLST not found')


def ST_to_genogroup(request,df):
  flag = False
  genogroup = str()
  associatedST = []
  for index,row in df.iterrows():
    if str(row['MLST']) == str(request):
      genogroup = row['Genogroup']
      flag = True
      break
  print(f'------------------------------------------------ \n')

  for index,row in df.iterrows():
    if row['Genogroup'] == genogroup:
      associatedST.append(row['MLST'])
  associatedST = sorted(associatedST)

  if flag == False:
    print(f'MLST {request} not found \n')
  else:
    print('Requested MLST found \n')
    if genogroup == 'singleST':
      print(f'MLST {request} belongs to singleST association')
      print(f'Genogroup cannot be determined')
      print(f'Viewing other single MLSTs: \n')
      output = str()
      count = 0
      for i in range(len(associatedST)):
        output += f'{associatedST[i]} '
        count += 1
        if len(output.split(' ')) == 11:
          print(output)
          output = str()
        elif str(output.split(' ')[-2]) == str(associatedST[-1]):
          print(output)
      print(f'{count} MLSTs total')

    elif len(associatedST) == 2:
      print(f'MLST {request} belongs to doubleST association')
      print(f'Genogroup cannot be determined')
      print(f'Viewing MLSTs with identical nucleotide profiles:\n ')
      for i in associatedST:
        print(i)
    else:
      print(f'MLST {request} belongs to genogroup {genogroup}')
      print(f'Viewing other MLSTs inside genogroup {genogroup}\n')
      output = str()
      count = 0
      for i in range(len(associatedST)):
        output += f'{associatedST[i]} '
        count += 1
        if len(output.split(' ')) == 4:
          print(output)
          output = str()
        elif str(output.split(' ')[-2]) == str(associatedST[-1]):
          print(output)
      print(f'{count} MLSTs total')
    print(f'------------------------------------------------ \n')





print('Before starting the program execution ensure that you have following packages installed and updated: Pandas, regex, biopython')
print("You can install them by terminating the program execution via EXIT and typing 'pip install -r requirements.txt\n")
print('-------------------------------------------------------------------------------------------------------------------------------')

while True:

  regime = str(input('Choose program mode: \n1.Obtain MLST from nucleotide profile\n2.Obtain nucleotide profile from MLST \n3.Determine which genogroup the sequence type belongs to\n4.Help and instruction\nEnter the required number to start program execution or get instruction \nEnter EXIT to cancel the program\nMode: '))


  if regime == '1':
    while True:
      request = str(input('Enter your nucleotide profile (capital ATGC letters only and N for doubtful positions; 18 SNPs). Enter RETURN to change the program mode: \n'))
      if request.upper() in ['RETURN', 'R']:
        break
      elif validation(request,regime) == False:
        print('Invalid request. Please, follow the instruction')
      elif validation(request,regime) == 'st_from_profile':
        Profile_to_ST (request,df)



  elif regime == '2':
    while True:
      request = str(input('Enter MLST. Enter RETURN to change the program mode: \n'))
      if request.upper() in ['RETURN', 'R']:
        break
      elif validation(request,regime) == False:
        print('Invalid request. Please, follow the instruction')
      elif validation(request,regime) == 'profile_from_st':
        ST_to_Profile(request,df)

  elif regime == '3':
    while True:
      request = str(input('Enter MLST. Enter RETURN to change the program mode: \n'))
      if request.upper() in ['RETURN', 'R']:
        break
      elif validation(request,regime) == False:
        print('Invalid request. Please, follow the instruction')
      elif validation(request,regime) == 'genogroup_from_st':
        ST_to_genogroup(request,df)

  elif regime.upper() in ['EXIT', 'E']:
    print('Program execution canceled')
    break

  elif regime == '4' or regime.upper() == 'HELP':
    help = str(f"------------------------------------------------------------------------------------------------------------------------------------------------------\n"
    f"Before starting the program execution ensure that you have following packages installed and updated: Pandas, regex, biopython\n"
    f"You can install them by terminating the program execution via EXIT and typing 'pip install -r requirements.txt\n"
    f"Database 'MLST_db_v2.csv' should be located in the current working directory before the program execution\n"
    f"For program navigation commands 'exit'/'e' and 'return'/'r' are allowed \n"
    f"Command 'exit' can be used for program termination\n"
    f"Command 'return' can be used to return to main menu\n"
    f"After program mode is chosen, we recommend to strictly follow terminal guidelines and insert appropriate data only to make the program execute correctly\n"
    f"For Linux console: "
    f"each of 3 program modes supports batch input if every element of the input array belongs to the new line (for example - past column from spreadsheet)\n"
    f"------------------------------------------------------------------------------------------------------------------------------------------------------\n" )
    print(help)


  else:
    print('Invalid regime. Please, follow the instruction')

