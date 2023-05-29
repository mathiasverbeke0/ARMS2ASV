import requests, re
import pandas as pd
import argparse

##################################
## PARSE COMMAND LINE ARGUMENTS ##
##################################

# Create parser
parser = argparse.ArgumentParser(description='Merging and subsequent MSA of multiple multifasta files')

# Add arguments
parser.add_argument('-e', '--ENA_accessions', required = True, help = 'The ENA accessions you want information for.')
parser.add_argument('-o', '--output', required = True, help = 'The output file.')

# Parse the arguments
args = parser.parse_args()

def get_coordinates(accession):
    # Construct the API URL
    api_url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
    
    # Make a request to the API
    response = requests.get(api_url)

    if response.status_code == 200:
        # Parse the XML response
        xml_data = response.text
        
        # Define the regular expression pattern
        pattern = "<DB>ENA-SAMPLE</DB>\s+<ID>(.*?)</ID>"

        # Find the ERS number using regular expression
        match = re.search(pattern, xml_data)

        # Extract the ERS number if a match is found
        if match:
            ers_number = match.group(1)
            # print("ERS number:", ers_number)

            # Construct the API URL
            api_url2 = f"https://www.ebi.ac.uk/ena/browser/api/xml/{ers_number}"
            
            # Make a request to the API
            response2 = requests.get(api_url2)

            if response.status_code == 200:
                # Parse the XML response
                xml_data2 = response2.text
                
                # Define the regular expression patterns
                pattern2 = "<TAG>lat lon</TAG>\s+<VALUE>(.*?), (.*?)</VALUE>"
                pattern3 = "<TAG>geographic location \(region and locality\)</TAG>\s+<VALUE>(.*?)</VALUE>"
                pattern4 = "<TAG>collection date</TAG>\s+<VALUE>(.*?)</VALUE>"
                pattern5 = "<TITLE>(.*?)</TITLE>"

                # Find the lon lat using regular expression
                match2 = re.search(pattern2, xml_data2)
                match3 = re.search(pattern3, xml_data2)
                match4 = re.search(pattern4, xml_data2)
                match5 = re.search(pattern5, xml_data2)

                # Extract the lon lat if a match is found
                if match2:
                    lat = match2.group(1)
                    lon = match2.group(2)

                else:
                    print("ERS number not found in the XML.")
                    lat = None
                    lon = None

                if match3:
                    region = match3.group(1)

                else:
                    print("Location not found in the XML.")
                    region = None

                if match4:
                    collection_date = match4.group(1)
                
                else:
                    print("Collection date not found in the XML.")
                    collection_date = None

                if match5:
                    title = match5.group(1)
                
                else:
                    print('Title not found in the XML.')
                    title = None

                return(lat, lon, region, collection_date, title)

            else:
                # Handle API request error
                print(f"API request failed with status code: {response.status_code}")
                return None, None, None, None, None

        else:
            print("ERS number not found in the XML.")
    
    else:
        # Handle API request error
        print(f"API request failed with status code: {response.status_code}")
        return None, None, None, None, None

accessions = args.ENA_accessions.split(',')

data = []

for accession in accessions:
    lat, lon, region, collection_date, title = get_coordinates(accession)
    data.append([lat, lon, region, collection_date, title])

df = pd.DataFrame(data, columns=['lat', 'long', 'region', 'collection_date', 'title'], index=accessions)

# Write the DataFrame to an Excel file

df.to_excel(args.output)
