import requests, string, time, sys
from bs4 import BeautifulSoup
import string
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys
import time

def get_species(prefix, incomplete):

    rSkips = 0
    species_names = []
    trials = 0

    while True: 
        
        if trials == 3:
            raise ValueError(f'Failed to establish a connection with WORMS after 3 attempts')
        
        url = f'http://www.marinespecies.org/aphia.php?p=taxlist&searchpar=0&tComp=begins&tName={prefix}&rSkips={rSkips}&vOnly=1&marine=1&fossil=4'
        response = requests.get(url)

        if response.status_code != 200:
            time.sleep(5)
            trials += 1
            continue

        # Parse the HTML
        soup = BeautifulSoup(response.content, 'html.parser')

        # Find the element containing the species information
        species_element = soup.find('ul', class_='list-group aphia_core_list_group_hover')

        if species_element == None and len(species_names) == 0:
            return [], incomplete
        
        elif species_element == None and len(species_names) != 0:
            return species_names, incomplete
        
        elif species_element != None and rSkips > 4900:
            incomplete.append(prefix)
            return species_names, incomplete

        else:
            # Extract the species names
            species_placeholder = [item.get_text(strip=True) for item in species_element.find_all('i')]
            species_placeholder = [item for item in species_placeholder if len(item.split()) == 2]
            
            species_names.extend(species_placeholder)
            rSkips += 100
            trials = 0

# Generate all possible prefixes from 'aaa' to 'zzz'
prefixes = [f"{c1}{c2}{c3}" for c1 in string.ascii_lowercase for c2 in string.ascii_lowercase for c3 in string.ascii_lowercase]

prefixes.extend(prefixes2)

# Iterate over each prefix and retrieve species names
all_species_names = []
incomplete = []

with open('/home/guest/Traineeship/ARMSProject/scripts/extra/extraspecies.txt', 'w') as file:

    # Use multithreading to send multiple requests to WORMS at the same time 
    with ThreadPoolExecutor(max_workers=16) as executor:
        # Submit tasks to the executor and store the future objects in a list
            futures = [executor.submit(get_species, prefix, incomplete) for prefix in prefixes]

            # Wait for all the tasks to complete and print the results
            for future in tqdm(as_completed(futures), total = len(prefixes)):
                try:
                    species_names, incomplete = future.result()
                    all_species_names.extend(species_names)
                    
                except Exception as e:
                    
                    print(f"An error occurred: {e}\nHalting all subsequent executions. This might take some time.")

                    # Cancel all remaining futures
                    for remaining_future in futures:
                        if not remaining_future.done():
                            remaining_future.cancel()
                    
                    sys.exit()

    # Convert the list to a set to remove duplicates
    all_species_names = list(set(all_species_names))
    
    print(f'{len(all_species_names)} unique species extracted from WORMS.')
    print(f'Prefixes with too many species to extract: {incomplete}')
    file.writelines(name + '\n' for name in all_species_names)