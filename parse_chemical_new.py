import re
import pandas as pd

class chemical_elements():
    def __init__(self, element):
        self.name = element
        self.A_site_r = None
        self.B_site_r = None
        self.Atomic_weight = None
        self.Volume_of_atom = None
        self.Nuclear_charge_effective_Slater = None
        self.Distance_from_core_electron_Schubert = None
        self.Distance_from_valence_electron_Schubert = None
        self.Electron_affinity = None
        self.Moment_nuclear_magnetic = None
        self.Valence_electron_number = None
        self.Pauling = None
        self.Allred_Rochow = None
        self._covalent_Martynov_Batsanov = None
        self.Absolute = None
        self.Oganov = None
        self.parameters = []

    def write_parameters(self, parameters: list):
        if len(parameters) != 15:
            raise ValueError(f"Неверное количество параметров {parameters} {self.name}")
        self.A_site_r = parameters[0]
        self.B_site_r = parameters[1]
        self.Atomic_weight = parameters[2]
        self.Volume_of_atom = parameters[3]
        self.Nuclear_charge_effective_Slater = parameters[4]
        self.Distance_from_core_electron_Schubert = parameters[5]
        self.Distance_from_valence_electron_Schubert = parameters[6]
        self.Electron_affinity = parameters[7]
        self.Moment_nuclear_magnetic = parameters[8]
        self.Valence_electron_number = parameters[9]
        self.Pauling = parameters[10]
        self.Allred_Rochow = parameters[11]
        self._covalent_Martynov_Batsanov = parameters[12]
        self.Absolute = parameters[13]
        self.Oganov = parameters[14]
        self.parameters = [float(param) for param in parameters]

def find_simple_element_concentrations(compaund, concentration_compaund=1.0):
    '''
    Функция для поиска элементов и их концентраций в строке. Разбирает только примитивные выражения, Которые записываются в формате [Элемент][Число][Элемент2][Число2] и т.д. В строке не должно быть одинаковых элементов'''
    pattern = r'([A-Z][a-z]*)(\d*\.?\d+)?'
    matches = re.findall(pattern, compaund)
    
    result = {}
    for element, concentration in matches:
        if not element:
            continue
        if concentration:
            conc_value = float(concentration.replace(',', '.'))*concentration_compaund
        else:
            conc_value = 1.0*concentration_compaund
        result[element] = conc_value
    
    return result

def get_compounds(formula):
    temp_formula = re.sub(
        r'\((.*?)\)',
        lambda m: '(' + re.sub(r'[-−]', '\x00', m.group(1)) + ')',
        formula,
        flags=re.DOTALL
    )
    compounds = temp_formula.split('-')
    compounds = [comp.replace('\x00', '-') for comp in compounds]
    return compounds

def get_sub_compaunds(compound: str, concentration_compaund=1.0) -> list:
    '''
    Находит сабкомпаунды в химической формуле.
    Сабкомпаундом называется вещество внутри скобок в формуле компаунда, 
    сабкомпаунд имеет свою концентрацию.
    
    Например:
    0.05Ca0.1(Bi0.5Na0.5)0.9ZrO3 вернет [['Bi0.5Na0.5', '0.9']]
    Если сабкомпаунда нет, то вернется исходное вещество
    
    Args:
        compound (str): Химическая формула
        
    Returns:
        list: Список списков [сабкомпаунд, концентрация] или False
    '''
    pattern = r'\(([^()]+)\)(\d*\.?\d*)'
    matches = re.findall(pattern, compound)
    
    if matches:
        result = []
        for sub_compound, concentration in matches:
            # Если концентрация не указана, считаем её равной 1
            concentration = concentration if concentration else '1'
            result.append([sub_compound, float(concentration)*concentration_compaund])
        return result
    else:
        return False
    
def rewrite_sub_compaunds(compound: str) -> str:
    '''
    Находит сабкомпаунды в химической формуле c его концентрацией, убирает скобки, изменя концентрацию входящих в него веществ.
    Сабкомпаундом называется вещество внутри скобок в формуле компаунда, 
    сабкомпаунд имеет свою концентрацию.
    
    Например:
    0.05Ca0.1(Bi0.5Na0.5)0.9ZrO3 - исходное вещество
    (Bi0.5Na0.5) - сабкомпаунд, концентрация 0.9
    вернет 0.05Ca0.1Bi0.45Na0.45ZrO3
    
    Args:
        compound (str): Химическая формула
        
    Returns:
        str: Новая химическая формула
    '''
    # Паттерн для поиска сабкомпаундов с их концентрацией
    #pattern = r'\(([^)]+)\)(\d*\.?\d+)?'

    pattern = r'\(([^()]+)\)(\d*\.?\d*)'
    #matches = re.findall(pattern, compound)
    
    while True:
        match = re.search(pattern, compound)
        if not match:
            break
            
        sub_compound = match.group(1)
        sub_concentration = float(match.group(2)) if match.group(2) else 1.0
        
        # Разбираем элементы в сабкомпаунде с учетом их концентрации
        elements = find_simple_element_concentrations(sub_compound, sub_concentration)
        
        # Собираем новую строку для замены
        replacement = ''
        for element, conc in elements.items():
            # Форматируем концентрацию без незначащих нулей
            formatted_conc = str(conc).rstrip('0').rstrip('.') if '.' in str(conc) else str(conc)
            replacement += f'{element}{formatted_conc}'
            
        # Заменяем сабкомпаунд в исходной строке
        compound = compound[:match.start()] + replacement + compound[match.end():]
    
    return compound

    


def reformat_formula(formula):
    formula = formula.replace(' ', '')
    formula = formula.replace('−', '-')
    formula = formula.replace('–', '-')
    formula = formula.replace(',', '.')
    return formula

def get_all_compaunds(formula, out_concentration=1.0):
    compaunds = get_compounds(formula)
    result= {}
    for compaund in compaunds:
        concentration1, compaund = get_concentrations_compaunds(compaund)
        inner_compaunds = get_compounds(compaund)
        if inner_compaunds != [compaund]:
            result={**result, **get_all_compaunds(compaund, concentration1*out_concentration)}
        else:
            result[compaund] = concentration1*out_concentration

    return result

def remove_outer_brackets(s: str) -> str:
    if len(s) < 2:
        return s
    
    if s[0] == '(' and s[-1] == ')':
        # Проверяем, что скобки действительно внешние (сбалансированы)
        balance = 0
        for i, char in enumerate(s):
            if char == '(':
                balance += 1
            elif char == ')':
                balance -= 1
                if balance == 0 and i != len(s) - 1:
                    # Нашли закрывающую скобку, которая не является последней
                    return s
        if balance == 0:
            return s[1:-1]
    
    return s

def get_concentrations_compaunds(compaund):
    match = re.match(r'^(\d+([.,]\d+)?)(.*)', compaund)
    if match:
        concentration = float(match.group(1).replace(',', '.'))
        formula = match.group(3)
        formula = remove_outer_brackets(formula)
        return [concentration, formula]
    else:
        # Если числа нет - концентрация 1, формула без изменений
        compaund = remove_outer_brackets(compaund)
        return [1.0, compaund]
def convert_to_dataframe(result):
    # Определяем максимальное количество пар элементов в A_site и B_site
    max_a = max((len(v["A_site"]))//2 for v in result.values()) if result else 0
    max_b = max((len(v["B_site"]))//2 for v in result.values()) if result else 0
    
    # Генерируем названия колонок
    columns = []
    for i in range(1, max_a + 1):
        columns.extend([f"A{i}", f"xA{i}"])
    for i in range(1, max_b + 1):
        columns.extend([f"B{i}", f"xB{i}"])
    
    # Собираем данные для DataFrame
    rows = []
    for formula, sites in result.items():
        row = [formula]
        
        # Обрабатываем A_site
        a_pairs = [(sites["A_site"][i], sites["A_site"][i+1]) 
                 for i in range(0, len(sites["A_site"]), 2)]
        for elem, x in a_pairs:
            row.extend([elem, x])
        # Заполняем недостающие элементы
        row.extend(["-", "-"] * (max_a - len(a_pairs)))
        
        # Обрабатываем B_site
        b_pairs = [(sites["B_site"][i], sites["B_site"][i+1]) 
                 for i in range(0, len(sites["B_site"]), 2)]
        for elem, x in b_pairs:
            row.extend([elem, x])
        # Заполняем недостающие элементы
        row.extend(["-", "-"] * (max_b - len(b_pairs)))
        
        rows.append(row)
    
    # Создаем DataFrame
    df = pd.DataFrame(rows, columns=["formula"] + columns)
    return df
def read_file(filename):
    with open(filename, 'r') as file:
        df = pd.read_excel(filename, engine='openpyxl')
        df = df.dropna(axis=1, how='all')
        return df
    
def generate_file_without_parameters(filename, dict_elements = None, test_element = None):
    if dict_elements == None:
        dict_elements = {}
    formulas = read_file(filename)
    result = {}
    len_b = 0
    len_a = 0
    count = 0
    dict_culculated_parameters = {}

                
    for col in formulas.columns:
        for formula in formulas[col]:
            is_debug = False
            if formula == test_element and test_element:
                print(f"исходная: {formula}")
                is_debug = True

            count += 1
            add_key = ""
            buf_formula = reformat_formula(formula)
            compaunds = get_all_compaunds(buf_formula)
            
            if is_debug:
                print("=======================================")
                print(f"{compaunds=}")
            if formula in result.keys():
                add_key = "(1)"
            else:
                add_key = ""
            #result[formula] = compaunds
            result[formula+add_key] = {
                "A_site": [],
                "B_site": []
            }
            dict_culculated_parameters[formula+add_key] = {
                "A_site": [],
                "B_site": []
            }
            if abs(1 - sum(compaunds.values())) > 0.07:
                result[formula+add_key]["A_site"].append(f"ошибка в стехиметрии составов {compaunds=} sum = {sum(compaunds.values())}")
                if is_debug:
                    print(f"{formula=} {compaunds=} sum = {sum(compaunds.values())}")
            else:
                for compaund, concentration in compaunds.items():
                    sub_compaund = get_sub_compaunds(compaund)
                    if is_debug:
                        print(f"{sub_compaund=}")

                    if sub_compaund:
                        compaund = rewrite_sub_compaunds(compaund)
                        if is_debug:
                            print(f"{compaund=}")

                    elements =  find_simple_element_concentrations(compaund)
                    if is_debug:
                        print(f"{elements=}")
                    A_site, B_site = split_perovskite(elements, concentration)
                    if is_debug:
                        print(f"{A_site=}")
                        print(f"{B_site=}")
                    for element, concentration in A_site.items():
                        result[formula+add_key]["A_site"].append(element)
                        result[formula+add_key]["A_site"].append(concentration)
                    for element, concentration in B_site.items():
                        result[formula+add_key]["B_site"].append(element)
                        result[formula+add_key]["B_site"].append(concentration)

            buf_parameters = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            for idx, element in enumerate(result[formula+add_key]["A_site"]):
                if element in dict_elements:
                    conc = float(result[formula+add_key]["A_site"][idx+1])
                    for index,param_val in enumerate(dict_elements[element].parameters):
                        buf_parameters[index] += conc*param_val
                else:
                    if isinstance(element, str) and dict_elements:
                        print(f"Для Элемента {element} не найдены параметры {formula=}")
                        
            dict_culculated_parameters[formula+add_key]["A_site"] = [round(param, 4) for param in buf_parameters]

            buf_parameters = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            for idx, element in enumerate(result[formula+add_key]["B_site"]):
                if element in dict_elements:
                    conc = float(result[formula+add_key]["B_site"][idx+1])
                    for index,param_val in enumerate(dict_elements[element].parameters):
                        buf_parameters[index] += conc*param_val
            dict_culculated_parameters[formula+add_key]["B_site"] = [round(param, 4) for param in buf_parameters]
            
    df = convert_to_dataframe(result)
    df.to_excel('result.xlsx')

    df = convert_results_to_dataframe(dict_culculated_parameters)
    df.to_excel('result_parameters.xlsx')

def convert_results_to_dataframe(result_dict):
    # Создаем заголовки столбцов
    param_names = [
    "A_site_r",
    "B_site_r",
    "Atomic_weight",
    "Volume_of_atom",
    "Nuclear_charge_effective(Slater)",
    "Distance_from_core_electron(Schubert)",
    "Distance_from_valence_electron(Schubert)",
    "Electron_affinity",
    "Moment_nuclear_magnetic",
    "Valence_electron_number",
    "Pauling",
    "Allred Rochow",
    "%_covalent_Martynov&Batsanov",
    "Absolute",
    "Oganov"
]
    a_columns = [f"{name}_A" for name in param_names]
    b_columns = [f"{name}_B" for name in param_names]
    columns = ['formula'] + a_columns + b_columns
    
    # Собираем данные
    data = []
    for formula, sites_data in result_dict.items():
        a_params = sites_data['A_site']
        b_params = sites_data['B_site']
        
        if len(a_params) != 15 or len(b_params) != 15:
            raise ValueError("Each site must contain exactly 15 parameters")
        
        row = [formula] + a_params + b_params
        data.append(row)
    
    return pd.DataFrame(data, columns=columns)

def split_perovskite(elements: dict, concentration_compaund=1.0, is_debug=False):
    if 'O' in elements:
        import copy
        elements2 = copy.deepcopy(elements)

        A_site_1 = {}
        B_site_1 = {}
        current_site1 = A_site_1
        summa_1 = 0
        for element, concentration in elements.items():
            if element == 'O':
                break
            current_site1[element] = concentration*concentration_compaund
            summa_1+= concentration
            if summa_1 >= 1:
                current_site1 = B_site_1
                summa_1 = 0


        A_site_2 = {}
        B_site_2 = {}
        current_site2 = B_site_2
        summa_2 = 0
        for element2, concentration2 in reversed( elements2.items() ):
            if element2 == 'O':
                continue

            current_site2[element2] = concentration2*concentration_compaund
            summa_2+= concentration2
            if summa_2 >= 1:
                current_site2 = A_site_2
                summa_2 = 0

        A_site_2 = dict(reversed(A_site_2.items()))
        B_site_2 = dict(reversed(B_site_2.items()))

        if is_debug:
            print(f"A_site_1: {A_site_1}")
            print(f"B_site_1: {B_site_1}")
            print(f"summa_1: {summa_1}")
            print(f"A_site_2: {A_site_2}")
            print(f"B_site_2: {B_site_2}")
            print(f"summa_2: {summa_2}")

        if abs(summa_1 - 1) > abs(summa_2 - 1):
            return A_site_2, B_site_2
        else:
            return A_site_1, B_site_1

    else:
        print(f"это не перовскит {elements} {concentration_compaund}")
        return {"это не перовскит":""}, {}



def test():
    result = {}
    formula = '(K0.44Na0.52Li0.04)0,995Cu0.0025(Nb0.86Ta0.1Sb0.04)O3'
    formula = reformat_formula(formula)
    compaunds = get_all_compaunds(formula)

    for compaund, concentration in compaunds.items():
        sub_compaund = get_sub_compaunds(compaund, concentration)
        if sub_compaund:
            compaund = rewrite_sub_compaunds(compaund)

        elements =  find_simple_element_concentrations(compaund)
        A_site, B_site = split_perovskite(elements, concentration, is_debug=True)
        #return A_site, B_site

def generate_file_with_parameters(filename, aboutelementsfile):

    elements_raw = read_file(aboutelementsfile)
    elements = {}
    #создали словарь с классами каждого элемента и заполненными аттрибутами
    for index,element in enumerate(elements_raw["element"]):
        elements[element] = chemical_elements(element)
        parameters = []
        for name in elements_raw.columns:
            if name == "element":
                continue
            parameters.append(elements_raw[name][index])
        elements[element].write_parameters(parameters)

    generate_file_without_parameters(filename, elements)

if __name__ == '__main__':
    filename = 'KNN_dataset_input.xlsx'
    aboutelementsfile = "elements_information.xlsx"

    #generate_file_with_parameters(filename, aboutelementsfile)
    generate_file_without_parameters(filename)
    #test()
    
            
