def switch(metal, bonded, distance=0):
    if 30 >= metal >= 21 and bonded == 8:
        if distance <= 1.7:
            return True
        return False
    elif 48 >= metal >= 39 and bonded == 8:
        if distance <= 1.91:    # TODO Fix this number - check Pd
            return True
        return False
    elif 80 >= metal >= 57 and bonded == 8:
        if distance <= 1.96:
            return True
        return False
    elif 30 >= metal >= 21 and bonded == 7:
        print("Row one to nitrogen")
        return False
    elif 48 >= metal >= 39 and bonded == 7:
        print("Row two bonded to nitrogen")
        return False
    elif 80 >= metal >= 57 and bonded == 7:
        print("Row three bonded to nitrogen")
        return False
    elif 30 >= metal >= 21 and bonded == 6:
        print("Row one to carbon")
        return False
    elif 48 >= metal >= 39 and bonded == 6:
        print("Row two bonded to carbon")
        return False
    elif 80 >= metal >= 57 and bonded == 6:
        print("Row three bonded to carbon")
        return False
