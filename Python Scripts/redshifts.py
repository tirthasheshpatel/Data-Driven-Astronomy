import numpy as np
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split

def get_features_targets(data):
    targets = data['redshift']
    features = np.zeros((data.shape[0], 4), dtype=np.float32)
    features[:,0] = data['u'] - data['g']
    features[:,1] = data['g']-data['r']
    features[:,2] = data['r']-data['i']
    features[:,3] = data['i']-data['z']
    return features, targets


if __name__ == "__main__":
    # load the data
    data = np.load('C:/Users/tirth/Desktop/Data Driven Astronomy/Python Scripts/sdss_galaxy_colors.npy')
    
    # call our function 
    features, targets = get_features_targets(data)
    features_train, features_test, targets_train, targets_test = train_test_split(features, targets)
    
    # print the shape of the returned arrays
    print(features_train[:2])
    print(targets_train[:2])

    model = DecisionTreeRegressor()
    model.fit(features_train, targets_train)
    preds = model.predict(features_test)
    print(1./targets.size * np.sum((targets_test-preds)**2))

    model_2 = SVR()
    model_2.fit(features_train, targets_train)
    preds_2 = model_2.predict(features_test)
    print(1./targets.size * np.sum((targets_test-preds_2)**2))