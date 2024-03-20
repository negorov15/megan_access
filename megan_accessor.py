# Main function of the package. Accesses MEGAN server and saves the datasets in csv file.
# Used packages: requests, HTTPBasicAuth
import os
import requests
from requests.auth import HTTPBasicAuth
# Import functions from other modules within package
from modificator import dataset_modifier, merge_data

# Function that accesses the megan server. Gets list of all data from the server
# Requires login and password and redirect the request to the function get_data(data).
# If status-code 200: proceeds with data
# If status-code 401: access denided. Invalid login or password.
# If status-code is other than 200 and 401: Error while accessing the server
def get_list():

    # Variables to store login and password
    login = input("login:")
    password = input("password:")
    # Variable saves respond from the server
    response = requests.get(url='http://maira.cs.uni-tuebingen.de:8001/megan6server/list',
                            auth=HTTPBasicAuth(login, password))
    # Variable stores the status code from server
    status_code = response.status_code

    if status_code == 200:
        # Store the response text and list all datasets available for analysis.
        data = response.text
        print("Available data: "+"\n"+data)
        return get_data(data, login, password)
    elif status_code == 401:
        print("Your login or password is invalid! Please try again!")
        return get_list()
    else:
        return "Error occurred while accessing the server. Status code: " + str(status_code)

# Function that allows to choose the files for analysis.
# Because the amount of data required for analysis is not fixed, user have to input "exit" command to proceed.
def get_data(data, login, password):

    # Splits the given response text in list.
    data_array = data.split()
    # List that stores all chosen datasets
    list_of_datasets = []
    print("Please, choose files for analysis. Once you have made your choice, print 'exit' ")
    # Infinite loop for choosing the datasets. Stops when the user inputs "exit"
    while True:
        usr_input = input("Choose a dataset: ")
        if (usr_input=="exit"):
            break
        if (usr_input not in data_array):
            print("Invalid input: no dataset in the list!")
            continue
        list_of_datasets.append(usr_input)
    return get_profile(list_of_datasets,login,password)

# Function that gets profiles for chosen datasets.
# Calls functions file_modifier(data) and create_profile(data)
def get_profile(datalist, login, password):
    # List with available classifications.
    classification = ['taxonomy', 'SEED', 'EGGNOG', 'GTDB', 'EC', 'INTERPRO2GO']
    print('Available classifications: ' + str(classification))
    num_of_datasets = len(datalist)
    while True:
        usr_input = input("What is the classification?: ")
        if usr_input not in classification:
            print('Error: classification is not available. Please enter a valid classification.')
        else: break
    for i in range(num_of_datasets):
        response = requests.get(
            url='http://maira.cs.uni-tuebingen.de:8001/megan6server' + "/getClassificationBlock?file=" + datalist[i] + "&classification=" + usr_input,
            auth=HTTPBasicAuth(login, password))
        sample_name = os.path.basename(datalist[i])
        datalist[i] = dataset_modifier(response.text, sample_name)
        #datalist[i] = dataset_modifier(response.text, 'Sample ' + str(i+1))
    return merge_data(datalist)

print(get_list())