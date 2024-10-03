



import subprocess
import sys
import importlib
import importlib.metadata
import platform

def install_and_import(package, module=None, version=None):
    if module is None:
        module = package
    try:
        installed_version = importlib.metadata.version(package)
        if version and installed_version != version:
            print(f"Updating {package} to version {version}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", f"{package}=={version}"])
    except importlib.metadata.PackageNotFoundError:
        print(f"Installing {package}...")
        if version:
            subprocess.check_call([sys.executable, "-m", "pip", "install", f"{package}=={version}"])
        else:
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    finally:
        # Import the module into the global namespace
        globals()[module] = importlib.import_module(module)




# List of required packages and their corresponding module names
required_packages = [
    ('pandas', None, '2.2.3'),
    ('xlsxwriter', None, '3.2.0'),
    ('requests', None, '2.28.1'),
    ('selenium', None, '4.4.0'),
    ('webdriver-manager', 'webdriver_manager', '4.0.2'),
    ('biopython', 'Bio', '1.83'),
]


for package, module, version in required_packages:
    install_and_import(package, module, version)


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
from Bio.Seq import Seq



#CONSTANTS
HELP = "\n _________________________________________________\n              | HIV Typing - COMET |                    \n _________________________________________________\n"
RULES_1 = "Maximum size of each sequence (excluding fasta) 20000 characters.\nThe sequences must be in fasta format, you can compress the file with gzip, max 50MB (8Mo).\n"
RULES_2 = "This script automatically renames the headers of the sequences to the assigned subtype. If you wish to disable this use --disable_auto_rename Y\nIt is advisable to run the following command in terminal before running the script:\npip install -r requirements.txt\n"
EXAMPLE_WINDOWS = "\nExample of usage:\n[Windows]\npython typing_using_comet.py --hiv_type 1 --fasta_file unknown_seqs.fasta --output_directory C:\\Users\\andre\\results --output_type xlsx\n"
EXAMPLE_MAC = "[MAC]\npython3 typing_using_comet.py --hiv_type 2 --fasta_file unknown_seqs.fasta.gz --output_directory C:\\Users\\andre\\results --output_type csv --disable_auto_rename Y\n"
OPTIONS_EXPLAINED_1 = "Options:\n--hiv_type [int] : 1 for HIV1, 2 for HIV2"
OPTIONS_EXPLAINED_2 = "--fasta_file [directory] : Directory of the fasta file (can be fasta.gz)"
OPTIONS_EXPLAINED_3 = "--output_directory [directory] : Directory for the output"
OPTIONS_EXPLAINED_4 = "--output_type [string] : \"xlsx\" for excel, \"csv\" for csv"
OPTIONS_EXPLAINED_5 = "--disable_auto_rename [optional] : \"Y\" for yes, otherwise omit from the commandline call\n"
ARGUMENT_OPTIONS = ["--hiv_type", "--fasta_file", "--output_directory", "--output_type", "--disable_auto_rename"]



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


def main(hiv_type, output_type, fasta_file_path, download_directory, disable_auto_rename):

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
                print(f"File saved as Excel to: {excel_file_path}")
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

def main_main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["Help", "hiv_type=", "fasta_file=", "output_directory=", "output_type=", "disable_auto_rename="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    hiv_type = None
    fasta_file = None
    output_directory = None
    output_type = None
    disable_auto_rename = None

    used_args = {"--hiv_type": None,
                 "--fasta_file": None,
                 "--output_directory": None,
                 "--output_type": None,
                 "--disable_auto_rename": None}

    for o, a in opts:
        if o in ("-h", "--Help"):
            print(HELP)
            print(RULES_1)
            print(RULES_2)
            print(OPTIONS_EXPLAINED_1)
            print(OPTIONS_EXPLAINED_2)
            print(OPTIONS_EXPLAINED_3)
            print(OPTIONS_EXPLAINED_4)
            print(OPTIONS_EXPLAINED_5)
            print(EXAMPLE_WINDOWS)
            print(EXAMPLE_MAC)
            sys.exit()
        elif o == "--hiv_type":
            hiv_type = a
            used_args[o] = a
        elif o == "--fasta_file":
            fasta_file = a
            used_args[o] = a
        elif o == "--output_directory":
            output_directory = a
            used_args[o] = a
        elif o == "--output_type":
            output_type = a
            used_args[o] = a
        elif o == "--disable_auto_rename":
            output_type = a
            used_args[o] = a
        else:
            assert False, "unhandled option"

    res = True
    for key in used_args.keys():
        if used_args[key] is None and key != "--disable_auto_rename":
            print(f"{key} argument is missing.")
            res = False

    if res == True:
        main(hiv_type, output_type, fasta_file, output_directory, disable_auto_rename)
    else:
        return


if __name__ == "__main__":
    main_main()



fasta_file_path_ = r"C:\Users\andre\OneDrive - IHMT-NOVA\helping_ihmt\victor_pim\HIV_Tools\Comet_for_Subtyping\files_for_testing\HIV2_B_seqs.fasta"
download_directory_ = r"C:\Users\andre\OneDrive - IHMT-NOVA\helping_ihmt\victor_pim\HIV_Tools\Comet_for_Subtyping"

# main(hiv_type=2, output_type="xlsx", fasta_file_path=fasta_file_path_, download_directory=download_directory_, disable_auto_rename="N")
