# Main function of the package. Accesses MEGAN server and saves the datasets in csv file.
import os
import requests
from requests.auth import HTTPBasicAuth
# Import functions from other modules within package
from modificator import dataset_modifier, merge_data


def get_list():
    """
    Grants access to the MeganServer, provided the user possesses the login and password,
    both of which are set to ‘guest‘ in this scenario.
    Displays all datasets provided by the active MeganServer instance.
    """
    # Variables to store login and password
    login = input("login:")
    password = input("password:")
    # Variable saves respond from the server
    response = requests.get(
        url="http://maira.cs.uni-tuebingen.de:8001/megan6server/list",
        auth=HTTPBasicAuth(login, password),
    )
    # Variable stores the status code from server
    status_code = response.status_code

    if status_code == 200:
        # Store the response text and list all datasets available for analysis.
        data = response.text
        print("Available data: " + "\n" + data)
        return get_data(data, login, password)
    elif status_code == 401:
        print("Your login or password is invalid! Please try again!")
        return get_list()
    else:
        return "Error occurred while accessing the server. Status code: " + str(
            status_code
        )


def get_data(data, login, password):
    """
    Allows the user to choose datasets for analysis.
    User can choose as many datasets as they want, but they have to be from the list of available datasets.
    Once the user has made their choice, the function calls get_profile() function to get profiles for chosen datasets.
    """
    # Splits the given response text in list.
    data_array = data.split()
    # List that stores all chosen datasets
    list_of_datasets = []
    print(
        "Please, choose files for analysis. Once you have made your choice, print 'exit' "
    )
    # Infinite loop for choosing the datasets. Stops when the user inputs "exit"
    while True:
        usr_input = input("Choose a dataset: ")
        if usr_input == "exit":
            break
        if usr_input not in data_array:
            print("Invalid input: no dataset in the list!")
            continue
        list_of_datasets.append(usr_input)
    return get_profile(list_of_datasets, login, password)


def get_profile(datalist, login, password):
    """
    Allows the user to choose the classification for analysis.
    User can choose only one classification, but it has to be from the list of available classifications.
    Once the user has made their choice, the function calls dataset_modifier() function to modify the datasets and merge_data() function to merge the datasets into a single dataset.
    """
    # List with available classifications.
    classification = ["taxonomy", "SEED", "EGGNOG", "GTDB", "EC", "INTERPRO2GO"]
    print("Available classifications: " + str(classification))
    num_of_datasets = len(datalist)
    while True:
        usr_input = input("What is the classification?: ")
        if usr_input not in classification:
            print(
                "Error: classification is not available. Please enter a valid classification."
            )
        else:
            break
    for i in range(num_of_datasets):
        response = requests.get(
            url="http://maira.cs.uni-tuebingen.de:8001/megan6server"
            + "/getClassificationBlock?file="
            + datalist[i]
            + "&classification="
            + usr_input,
            auth=HTTPBasicAuth(login, password),
        )
        sample_name = os.path.basename(datalist[i])
        datalist[i] = dataset_modifier(response.text, sample_name)
        #datalist[i] = dataset_modifier(response.text, 'Sample ' + str(i+1))
    return merge_data(datalist)


print(get_list())