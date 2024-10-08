




# Now import the modules
import getopt
import os
import time
import glob
import traceback
import sys
import requests
from urllib.parse import urljoin
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
import pandas as pd
from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse


#CONSTANTS
# RULES_1 = "Maximum size of each sequence (excluding fasta) 20000 characters.\nThe sequences must be in fasta format, you can compress the file with gzip, max 50MB (8Mo).\n"
# RULES_2 = "This script automatically renames the headers of the sequences to the assigned subtype. If you wish to disable this use --disable_auto_rename Y\nIt is advisable to run the following command in terminal before running the script:\npip install -r requirements.txt\n"
# EXAMPLE_WINDOWS = "\nExample of usage:\n[Windows]\npython typing_using_comet.py --hiv_type 1 --fasta_file unknown_seqs.fasta --output_directory C:\\Users\\results --output_type xlsx --split_per_subtype B G\n"
# EXAMPLE_MAC = "[MAC]\npython3 typing_using_comet.py --hiv_type 2 --fasta_file unknown_seqs.fasta.gz --output_directory C:\\Users\\results --output_type csv --disable_auto_rename Y --split_per_subtype C\n"
# OPTIONS_EXPLAINED_1 = "[int] : 1 for HIV1, 2 for HIV2"
# OPTIONS_EXPLAINED_2 = "[directory] : Directory of the fasta file (can be fasta.gz)"
# OPTIONS_EXPLAINED_3 = "[directory] : Directory for the output"
# OPTIONS_EXPLAINED_4 = "[string] : \"xlsx\" for excel, \"csv\" for csv"
# OPTIONS_EXPLAINED_5 = "[optional] : \"Y\" for yes, otherwise omit from the commandline call"
# OPTIONS_EXPLAINED_6 = "[string] : Can be one type or multiple characters. Eg. --split_per_subtype G B\n"
# ARGUMENT_OPTIONS = ["--hiv_type", "--fasta_file", "--output_directory", "--output_type", "--disable_auto_rename", "--split_per_subtype"]
#

HELP = """
_________________________________________________
              | HIV Typing - COMET |             
_________________________________________________
Description:
  This script automatically renames the headers of the sequences to the assigned subtype.
  If you wish to disable this use --disable_auto_rename
  It is advisable to run the following command in terminal before running the script:
    pip install -r requirements.txt.

Example of usage:
[WINDOWS]
  python script.py --hiv_type 2 --fasta_file data.fasta --output_directory ./output --output_type xlsx --split_per_subtype A B C

[MAC]
  python script.py --hiv_type 2 --fasta_file data.fasta --output_directory ./output --output_type csv --disable_auto_rename Y --split_per_subtype G

Required Arguments:
  --hiv_type [int]                      Specify the HIV type (e.g., 1 for HIV1, 2 for HIV2).
  --fasta_file [directory]              Directory of the fasta file (can be fasta.gz).
  --output_directory [directory]        Directory to save output files.
  --output_type [string]                Specify the output file type (e.g., csv or xlsx). Default is xlsx

Optional Arguments:
  --disable_auto_rename [string]        Disable automatic renaming of output files, Y for yes, otherwise omit from the commandline call
  --split_per_subtype [string]          Single or multiple subtypes for splitting per subtype analysis, (e.g --split_by_subtype A or --split_by_subtype B C)

Rules:
  - Maximum size of each sequence (excluding fasta) 20000 characters.
  - The sequences must be in fasta format, you can compress the file with gzip, max 50MB (8Mo)

For more information, refer to the documentation.
"""


def wait_for_download_completion(download_path, timeout=60):
    seconds = 0
    while seconds < timeout:
        time.sleep(1)
        files = os.listdir(download_path)
        if any(filename.endswith('.crdownload') for filename in files):
            print("Waiting for download to complete...")
        else:
            return True
        seconds += 1
    return False


def get_latest_file(download_path):
    list_of_files = glob.glob(os.path.join(download_path, '*'))
    if not list_of_files:
        return None
    latest_file = max(list_of_files, key=os.path.getctime)
    return latest_file


def read_fasta(filename: str):
    records = list(SeqIO.parse(filename, "fasta"))
    return records


def main(hiv_type, output_type, fasta_file_path, download_directory, disable_auto_rename, wanted_subtypes):

    print("\n\n\tYour options:")
    # Your main function implementation goes here
    print("hiv_type:", hiv_type)
    print("output_type:", output_type)
    print("fasta_file:", fasta_file_path)
    print("output_directory:", download_directory)
    print("disable_auto_rename:", disable_auto_rename)
    print("subtypes:", wanted_subtypes)



    # Ensure the output directory exists
    global df
    if not os.path.exists(download_directory):
        os.makedirs(download_directory)

    # Create the output directory exists
    output_directory = os.path.join(download_directory, "Comet_Output")
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Extract the base filename without directory
    input_filename = os.path.basename(fasta_file_path)

    # Set up Chrome options
    chrome_options = Options()
    chrome_prefs = {
        'download.default_directory': download_directory,
        'download.prompt_for_download': False,
        'download.directory_upgrade': True,
        'safebrowsing.enabled': True,
        'profile.default_content_settings.popups': 0,
        'profile.default_content_setting_values.automatic_downloads': 1,
        'safebrowsing.disable_download_protection': True,  # Use with caution
    }
    chrome_options.add_experimental_option('prefs', chrome_prefs)

    # Uncomment to run in headless mode
    # chrome_options.add_argument('--headless')

    # Set up the WebDriver with options
    # driver = webdriver.Chrome(ChromeDriverManager().install()) # OLD version, works on WINDOWS, not on MAC

    from selenium.webdriver.chrome.service import Service

    # Initialize the WebDriver with service and options
    driver = webdriver.Chrome(
        service=Service(ChromeDriverManager().install()),
        options=chrome_options
    )

    print()

    try:
        if hiv_type == "1":
            # Navigate to the COMET HIV-1 page
            driver.get('https://comet.lih.lu/index.php?cat=hiv1')
            print(f"Navigated to {driver.current_url}")
        elif hiv_type == "2":
            # Navigate to the COMET HIV-1 page
            driver.get('https://comet.lih.lu/index.php?cat=hiv2')
            print(f"Navigated to {driver.current_url}")

        # Wait for the file input to be present
        file_input = WebDriverWait(driver, 60).until(
            EC.presence_of_element_located((By.NAME, 'fastafile'))
        )
        print("File input found.")


        file_input.send_keys(fasta_file_path)
        print(f"Uploaded file: {fasta_file_path}")

        # Check the non-commercial checkbox
        checkbox = driver.find_element(By.NAME, 'non_commercial')
        checkbox.click()
        print("Non-commercial checkbox checked.")

        # Click the submit button
        submit_button = driver.find_element(By.NAME, 'submit')
        submit_button.click()
        print("Submit button clicked.")

        # Define the expected link text
        expected_link_text = "Download the results in CSV format (tab separated, available for 3 days)"

        # Wait for the link with the exact text to be present and clickable
        download_link = WebDriverWait(driver, 120).until(
            EC.text_to_be_present_in_element(
                (By.XPATH, '//a[contains(@href, "/csv.php?job=")]'),
                expected_link_text
            )
        )
        print("Download link with expected text found.")

        # Now locate the link element
        download_link_element = driver.find_element(By.XPATH,
                                                    f'//a[contains(@href, "/csv.php?job=") and text()="{expected_link_text}"]')

        # Extract the href attribute to get the download URL
        download_href = download_link_element.get_attribute('href')
        print(f"Download href: {download_href}")

        # Construct the full download URL
        base_url = 'https://comet.lih.lu'
        download_url = urljoin(base_url, download_href)
        print(f"Download URL: {download_url}")


        # Option 2 (Alternative): Download the file using requests
        # Get cookies from Selenium and add them to requests session
        session = requests.Session()
        for cookie in driver.get_cookies():
            session.cookies.set(cookie['name'], cookie['value'])
        # Download the file
        response = session.get(download_url)


        # After downloading the file using requests
        if response.status_code == 200:
            if output_type == "xlsx":
                # Construct the output filename
                excel_file_path = os.path.join(output_directory, f'{os.path.splitext(input_filename)[0]}_results.xlsx')

                # Read the CSV content into a pandas DataFrame
                # Since the data is in CSV format, we'll use StringIO to read it from the response content
                csv_content = response.content.decode('utf-8')
                df = pd.read_csv(StringIO(csv_content), sep='\t')  # Adjust the separator if necessary

                df = df.rename(columns={"name": "sequence"})

                # Save the DataFrame as an Excel file
                df.to_excel(excel_file_path, index=False)
                print(f"File saved as Excel to: {excel_file_path}\n")
            elif output_type == "csv":
                # Construct the output filename
                csv_file_path = os.path.join(output_directory, f"{os.path.splitext(input_filename)[0]}_results.csv")
                # Save the file
                with open(csv_file_path, 'wb') as f:
                    f.write(response.content)
                print(f"File downloaded via requests to: {csv_file_path}")
        else:
            print(f"Failed to download file via requests. Status code: {response.status_code}")


    except Exception as e:
        print(f"An error occurred: {e}")
        traceback.print_exc()
        driver.save_screenshot('error_screenshot.png')
        with open('error_page_source.html', 'w', encoding='utf-8') as f:
            f.write(driver.page_source)
    finally:
        # Close the browser
        driver.quit()

    print("\nRenaming headers in Fasta")
    if disable_auto_rename != "Y":
        records = read_fasta(fasta_file_path)

        for record in records:
            if record.id in df["sequence"].to_list():
                found_type = df.loc[df["sequence"] == record.id, "subtype"].values[0]
                record.id = found_type + "." + record.id

        renamed_fasta_dir = os.path.join(output_directory, f"{os.path.splitext(input_filename)[0]}_typed.fasta")
        with open(renamed_fasta_dir, "w") as output_handle:

            SeqIO.write(records, output_handle, "fasta")
        print(f"Renamed headers, file found in {renamed_fasta_dir}\n")

    print(f"Filtering for {wanted_subtypes} subtypes\n")
    filename = f"{os.path.splitext(input_filename)[0]}_typed.fasta".split(".")[0]

    records = read_fasta(renamed_fasta_dir)

    all_record_ids = []

    filtered_records = []

    for record in records:
        for wanted_subtype in wanted_subtypes:
            if wanted_subtype == record.id.split(".")[0]:
                wanted_record = SeqRecord(
                    Seq(f"{record.seq}"),
                    id=record.id)
                filtered_records.append(wanted_record)
        all_record_ids.append(record.id.split(".")[0])

    not_in_identified = [wanted_subtype for wanted_subtype in wanted_subtypes if wanted_subtype not in all_record_ids]

    kept_subtypes = list(set(wanted_subtypes) - set(not_in_identified))

    if len(not_in_identified) == 0:
        print("All wanted subtypes are present at least once in the filtered file.\n")
    else:
        print(f"Subtypes {set(not_in_identified)} are not present in the filtered file.\n")

    write_result = SeqIO.write(filtered_records, os.path.join(output_directory, f"{filename}_filtered.fasta"), "fasta")
    print(f"Filtered fasta written with {kept_subtypes} subtypes.\nFile found in: ",  os.path.join(output_directory, f"{filename}_filtered.fasta\n"))


def main_main():
    # Check if help was requested before parsing arguments
    if '-h' in sys.argv or '--Help' in sys.argv:
        # Display your custom help messages
        print(HELP)
        sys.exit()

    parser = argparse.ArgumentParser(description='Process HIV sequence data.', add_help=False)

    # Required arguments
    parser.add_argument('--hiv_type', required=True, help='Specify the HIV type.')
    parser.add_argument('--fasta_file', required=True, help='Path to the input FASTA file.')
    parser.add_argument('--output_directory', required=True, help='Directory to save output files.')
    parser.add_argument('--output_type', required=True, help='Specify the output file type.')

    # Optional arguments
    parser.add_argument('--disable_auto_rename', action='store_true', help='Disable automatic renaming of output files.')
    parser.add_argument('--split_per_subtype', nargs='+', help='List of subtypes for splitting per subtype analysis.')

    # Custom help option
    parser.add_argument('-h', '--Help', action='store_true', help='Show this help message and exit.')

    args = parser.parse_args()

    # Assign arguments to variables
    hiv_type = args.hiv_type
    fasta_file = args.fasta_file
    output_directory = args.output_directory
    output_type = args.output_type
    disable_auto_rename = args.disable_auto_rename
    subtypes = args.split_per_subtype if args.split_per_subtype else []

    # Call the main function with the parsed arguments
    main(hiv_type, output_type, fasta_file, output_directory, disable_auto_rename, subtypes)

if __name__ == "__main__":
    main_main()