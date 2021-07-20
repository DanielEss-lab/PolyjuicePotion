import re


class ChargeChanger:
    def __init__(self, filename):
        self.filename = filename

    def change(self, charge):
        with open(self.filename, 'r') as file:
            filedata = file.read()

        old_charge = int(re.findall("q = (-?\d+)", filedata)[0])
        new_charge = old_charge + charge

        replace_string = "q = " + str(old_charge)
        new_string = "q = " + str(new_charge)
        filedata = filedata.replace(replace_string, new_string)

        with open(self.filename, 'w') as file:
            file.write(filedata)


if __name__ == "__main__":
    changer = ChargeChanger("W3095.xyz")
    changer.change(0)
