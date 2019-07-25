
import requests
import pandas as pd
import os

#%% User Input
username='eCommons_id'
password='MyPassword'


#%%
r = requests.post('https://reagenttracker.hms.harvard.edu/api/v0/login',
                  json={'username': username, 'password': password})
cookies = r.cookies
rt_response = requests.get('https://reagenttracker.hms.harvard.edu/api/v0/search?q=', cookies=cookies).json()
sm = [d for d in rt_response['canonicals'] if d['type'] == 'small_molecule' and d['name'] != 'DEPRECATED']
#%%
rt_data = pd.DataFrame(sm)[['lincs_id', 'name', 'alternate_names', 'smiles']]
#rt_data.set_index('lincs_id', inplace=True, verify_integrity=True)
#rt_data.rename(columns={'smiles': 'smiles_original'}, inplace=True)
#rt_data['smiles_original'].fillna('', inplace=True)
rt_data['smiles'].fillna('', inplace=True)


print(rt_data[0:5])
#%%
dir_input='/users/nienkemoret/Dropbox (HMS)/CBDM-SORGER/Collaborations/LSP_data_organization/ChemInformatics_Files/Reagenttracker'
os.chdir(dir_input)

molecule_list_df=pd.DataFrame(rt_data)
#molecule_list_df.columns=['hms_id','pref_name','alt_name','smiles']

molecule_list_df.to_csv('smallmol_list_RT_20190617.csv', index=False)