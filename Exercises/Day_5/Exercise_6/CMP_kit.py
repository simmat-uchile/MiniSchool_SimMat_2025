import numpy as np
from numpy.linalg import norm
from scipy import stats
from scipy import ndimage
from scipy.interpolate import CubicSpline
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import sys
import math
import os
from scipy.interpolate import interp1d

###funciones auxiliares. 

def read_file_info(FILE1):
    NATOMS=0 
    GRID, POINTS = None, None
    alattvec, blattvec, clattvec = None, None, None
    
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i == 2:
                alattvec = list(map(float, line.split()))
            elif i == 3:
                blattvec = list(map(float, line.split()))
            elif i == 4:
                clattvec = list(map(float, line.split()))
            elif i == 6:
                line = map(int, line.split())
                NATOMS = sum(line)
            elif i == NATOMS + 9:
                GRID = list(map(int, line.split()))
                POINTS = np.prod(GRID)
            elif i > NATOMS + 11:
                break
    SKIPELINES = NATOMS + 9
    if POINTS % 5 == 0:
        MAXROWS = int(POINTS/5)
    elif POINTS % 5 != 0:
        MAXROWS = int(np.ceil(POINTS/5))
    print('NATOMS = ', NATOMS, 'GRID = ',GRID, 'POINTS = ',POINTS)
    return NATOMS, GRID, POINTS, alattvec, blattvec, clattvec, SKIPELINES, MAXROWS

def process_file(filename, skip_lines, max_rows, points):
    df = pd.read_fwf(filename, colspecs='infer', skiprows=skip_lines, nrows=max_rows)
    data = df.to_numpy().flatten('C')[:points]
    return data

def calculate_fukui_and_r2(dn, CHG1, CHG2, CHG3, CHG4):
    density_matrix = np.vstack([CHG1, CHG2, CHG3, CHG4]).astype(float)
    coeffs = np.polyfit(dn, density_matrix, deg=1) 
    slopes = coeffs[0]  
    correlation_matrix = np.corrcoef(density_matrix)
    r_squared_values = correlation_matrix[0, 1] ** 2

    return slopes, r_squared_values

def write_fukui_file(FILE1, NATOMS, fukui, filename):
    ### Abrir el archivo de salida para escribir todo de una vez
    with open(filename, "w") as FUKUIFILE:
        
        # Escribir el encabezado desde FILE1
        with open(FILE1) as fp:
            for i, line in enumerate(fp):
                if i < NATOMS + 10:
                    FUKUIFILE.write(line)
                else:
                    break

        # Ahora procesamos y escribimos los datos de 'fukui'
        num_full_rows = fukui.size // 5
        last_row_size = fukui.size % 5
        
        if last_row_size == 0:
            fukui_matrix = fukui.reshape(-1, 5)
        else:
            fukui_matrix = np.zeros((num_full_rows + 1, 5))
            fukui_matrix[:num_full_rows] = fukui[:num_full_rows * 5].reshape(-1, 5)
            fukui_matrix[-1, :last_row_size] = fukui[num_full_rows * 5:]
        
        # Escribir los datos de 'fukui' en el archivo
        for row in fukui_matrix:
            formatted_row = " ".join(f"{value: .11E}" for value in row if value != 0)
            FUKUIFILE.write(formatted_row + "\n")
    print('file saved: ',filename)
    return filename

def reshape_xyz(NGX,NGY,NGZ, CHG):
    GHGtem =  np.zeros((NGX,NGY,NGZ))
    for N3 in range(NGX):
        for N2 in range(NGY):
            for N1 in range(NGZ):
                GHGtem[N3,N2,N1] = CHG[N1*NGX*NGY + N2*NGX + N3]
    return GHGtem

def missing(data, POINTS):
    nan_count=data.isna().sum().sum()  
    if nan_count==0:
        CHG=data.to_numpy()
        CHG=CHG.flatten('C')
        CHG=CHG[0:POINTS]
        return CHG

    else:
        CHG=data.to_numpy()
        CHG=CHG.flatten('C')
        CHG=CHG[0:POINTS]
        hay_nan = any(np.isnan(valor) for valor in CHG)
        if hay_nan:
            print("Input error")
        else:
            print("")
        return CHG

def func_correction(omega_vol, eps, LS, LZ, NGZ):
        q = 1.6022 *1e-19
        eps_0 = 8.854*1e-22 
        modulo = (-q/(2.0*eps_0*omega_vol)) 
        correction = []
        recorrido = np.linspace(-(LZ*0.5),(LZ*0.5) , NGZ)
        for z in recorrido:
            if (-LS*0.5) < z < (LS*0.5):
                func_dentro = (1/eps)*(z*z + LS*LS*0.25*(eps -1) - eps*LZ*LZ*0.25)
                correction.append(func_dentro*modulo)
            elif (abs(z)>= (LS*0.5)):
                func_fuera = z*z - (LZ*LZ*0.25)
                correction.append(func_fuera*modulo)
            else:
                print("Error func_correction") ### aqui hubo una i
        return correction

def vector_atoms_position(nombre_archivo, N_ATOMS):
        fin = 8 + N_ATOMS
        with open(nombre_archivo, 'r') as archivo:
            lineas = archivo.readlines()
            lineas_interesantes = lineas[8 : fin]
            vectores = [list(map(float, linea.strip().split())) for linea in lineas_interesantes]
            return vectores

def compute_lattice_parameters(alattvec, blattvec, clattvec):
    LATTCURA = (1.0 / 0.529177210903) * np.dstack([alattvec, blattvec, clattvec])
    LATTCURA = LATTCURA[0]
    LATTCURB = np.zeros((3, 3))
    LATTCURB[0] = np.cross(LATTCURA[1], LATTCURA[2])
    LATTCURB[1] = np.cross(LATTCURA[2], LATTCURA[0])
    LATTCURB[2] = np.cross(LATTCURA[0], LATTCURA[1])
    VOL = np.abs(np.linalg.det(np.dstack([alattvec, blattvec, clattvec])))
    omega = np.dot(np.cross(LATTCURA[1], LATTCURA[2]), LATTCURA[0])
    LATTCURB /= omega
    return LATTCURB, VOL, omega

def compute_gsquared(GRID, LATTCURB, q):
    NGX, NGY, NGZ = GRID
    LPCTX = np.array([NX if NX < int(NGX / 2) + 1 else NX - NGX for NX in range(NGX)])
    LPCTY = np.array([NY if NY < int(NGY / 2) + 1 else NY - NGY for NY in range(NGY)])
    LPCTZ = np.array([NZ if NZ < int(NGZ / 2) + 1 else NZ - NGZ for NZ in range(NGZ)])
    GSQU = np.zeros((NGX, NGY, NGZ))
    
    for N3 in range(NGZ):
        for N2 in range(NGY):
            for N1 in range(NGX):
                GX = LPCTX[N1] * LATTCURB[0, 0] + LPCTY[N2] * LATTCURB[0, 1] + LPCTZ[N3] * LATTCURB[0, 2]
                GY = LPCTX[N1] * LATTCURB[1, 0] + LPCTY[N2] * LATTCURB[1, 1] + LPCTZ[N3] * LATTCURB[1, 2]
                GZ = LPCTX[N1] * LATTCURB[2, 0] + LPCTY[N2] * LATTCURB[2, 1] + LPCTZ[N3] * LATTCURB[2, 2]
                GSQU[N1, N2, N3] = (GX * GX + GY * GY + GZ * GZ)
    
    GSQU[0, 0, 0] = 1.0
    GSQU = 1.0 / (GSQU * (2.0 * np.pi * 2.0 * np.pi))
    GSQU[0, 0, 0] = q
    return GSQU

def compute_planar_average_nz(CHG00, NGX, NGY, NGZ, axis):
    if axis == 'x':
        return [np.sum(CHG00[nx, :, :]) / (NGY * NGZ) for nx in range(NGX)]
    if axis == 'y':
        return [np.sum(CHG00[:, ny, :]) / (NGX * NGZ) for ny in range(NGY)]
    if axis == 'z':
        return [np.sum(CHG00[:, :, nz]) / (NGY * NGX) for nz in range(NGZ)]

def plot_planar_average(z_axis, PLANAR1, PLANAR3):
    plt.plot(z_axis, PLANAR3, label='Corrected', linewidth=2)
    plt.plot(z_axis, PLANAR1, label='No correction', linewidth=2)
    plt.title('Planar Average Fukui Potential', fontsize=18)
    plt.ylabel(r'$v_{f}(r)$ (eV)', fontsize=12.5)
    plt.xlabel(r'Z-direction ($\AA$)', fontsize=12.5)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.legend()
    plt.savefig('PA_vFuk.svg')

def read_data(FILE0, skipelines, maxrows, POINTS):
    if (POINTS % 5) != 0:
      with open(FILE0, 'r') as file:
         last_line = file.readlines()[maxrows + skipelines]
         df0 = pd.read_table(FILE0, sep=r'\s+', skiprows=skipelines+1, names=range(5), nrows=maxrows)
    else:
        skiprowsn=skipelines + 1
        df0 = pd.read_table(FILE0, sep=r'\s+', skiprows=skipelines+1, names=range(5), nrows=maxrows)
    return df0

def closest_value_position(lst, number, percentage=0.4):
    num_elements = int(len(lst) * percentage)
    position, closest_value = min(enumerate(lst[-num_elements:]), key=lambda x: abs(x[1] - number))
    return closest_value, position + len(lst) - num_elements

def extract_data_from_file(file_path, marker='#z-vcor.dat'):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Find the start of the data after the marker
    data_start = None
    for i in range(len(lines)-1, -1, -1):
        if lines[i].strip().startswith(marker):
            data_start = i + 1
            break

    if data_start is None:
        raise ValueError(f"The file does not contain the marker '{marker}'")

    # Extract the data starting from the marker
    data_lines = lines[data_start:]
    z2, value2 = [], []
    for line in data_lines:
        if line.strip():
            cols = line.split()
            z2.append(float(cols[0]))
            value2.append(float(cols[1]))

    # Convert to numpy arrays
    return np.array(z2), np.array(value2)

def write_xyz(CHG, NGX, NGY, NGZ, LATTCURADO, omega, output_file='xyz_value.dat'):
    suma = 0.0
    with open(output_file, 'w') as xyz_file:
        # Write the header for the XYZ file
        xyz_file.write(f"{NGX * NGY * NGZ}\n")
        xyz_file.write("Generated XYZ file\n")

        # Iterate over the grid and calculate coordinates and charge density
        for N3 in range(NGZ):
            for N2 in range(NGY):
                for N1 in range(NGX):
                    # Calculate the X, Y, Z coordinates
                    x = (N1 * LATTCURADO[0, 0] + N2 * LATTCURADO[1, 0] + N3 * LATTCURADO[2, 0]) / NGX
                    y = (N1 * LATTCURADO[0, 1] + N2 * LATTCURADO[1, 1] + N3 * LATTCURADO[2, 1]) / NGY
                    z = (N1 * LATTCURADO[0, 2] + N2 * LATTCURADO[1, 2] + N3 * LATTCURADO[2, 2]) / NGZ

                    # Get the charge density value at the current grid point
                    valor_densidad = CHG[N1, N2, N3]
                    
                    # Accumulate the sum of charge density values
                    suma += valor_densidad

                    # Write the coordinates and charge density to the file
                    xyz_file.write(f"{x:>12.4f}{y:>12.4f}{z:>12.4f}{valor_densidad:>20.8e}\n")

    # Calculate the differential element
    Elemento_volumen = omega / (NGX * NGY * NGZ)
    
    # Calculate the integral of the charge density
    Integral_rho = suma * Elemento_volumen
    return output_file,Integral_rho

def detect_local_extrema_3D(data, order, extrema_type='min'):
    """
    Detects local extrema (maxima or minima) in a 3D array.

    Parameters
    ----------
    data : ndarray
        3D array of data.
    order : int
        Number of points on each side to use for comparison.
    extrema_type : str
        Type of extrema to detect ('min' for minima, 'max' for maxima).

    Returns
    -------
    coords : ndarray
        Coordinates of the local extrema.
    values : ndarray
        Values of the local extrema.
    """
    size = 1 + 2 * order
    footprint = np.ones((size, size, size))
    footprint[order, order, order] = 0

    if extrema_type == 'min':
        # Detect local minima
        filtered = ndimage.minimum_filter(data, footprint=footprint, mode='wrap')
        mask_extrema = data < filtered
    elif extrema_type == 'max':
        # Detect local maxima
        filtered = ndimage.maximum_filter(data, footprint=footprint, mode='wrap')
        mask_extrema = data > filtered
    else:
        raise ValueError("extrema_type must be 'min' or 'max'")

    coords = np.asarray(np.where(mask_extrema)).T
    values = data[mask_extrema]

    return coords, values

def write_formatted_data_localm(filename, data):
    with open(filename, 'w') as archivo:
        for fila in data:
            archivo.write(f"{fila[0]:>10.7f} {fila[1]:>10.7f} {fila[2]:>10.7f} {fila[3]:>10.7f}\n")
    print(f"Data successfully written to {filename}")

def calcular_xyz_val(coords, values, GRID, alattvec, blattvec, clattvec):
    """
    Función para calcular la colección de min_xyz_val para una serie de coordenadas y valores.

    Args:
    - coords (list or np.ndarray): Lista o array de coordenadas [NZ, NY, NX].
    - values (list or np.ndarray): Lista o array de valores asociados a cada coordenada.
    - GRID (list or np.ndarray): Tamaño de la cuadrícula [GRID_X, GRID_Y, GRID_Z].
    - alattvec, blattvec, clattvec (list or np.ndarray): Vectores de red asociados.

    Returns:
    - min_xyz_vals (list): Lista de arrays que contienen [x_min, y_min, z_min, value] para cada coordenada.
    """
    
    min_xyz_vals = []

    for i in range(len(values)):
        NZ = coords[i][0]
        NY = coords[i][1]
        NX = coords[i][2]

        x_min = ((NX - 1) / GRID[0]) * np.linalg.norm(np.array(alattvec), 2)
        y_min = ((NY - 1) / GRID[1]) * np.linalg.norm(np.array(blattvec), 2)
        z_min = ((NZ - 1) / GRID[2]) * np.linalg.norm(np.array(clattvec), 2)

        min_xyz_val = np.array([x_min, y_min, z_min, values[i]])

        min_xyz_vals.append(min_xyz_val)

    return min_xyz_vals

def convert_locpot(LOCPOT_cor2, NGX, NGY, NGZ):
    LOCPOTtem = np.zeros(NGX * NGY * NGZ)
    for N1 in range(NGZ):
        for N2 in range(NGY):
            for N3 in range(NGX):
                LOCPOTtem[N1 * NGX * NGY + N2 * NGX + N3] = - LOCPOT_cor2[N3, N2, N1]
                #cambiamos el signo al tiro
    return LOCPOTtem

def read_xyzval(archivo):
    datos = []
    with open(archivo, 'r') as f:
        for linea in f:
            try:
                x, y, z, valor = map(float, linea.strip().split()[:4])
                datos.append((x, y, z, valor))
            except ValueError:
                continue
    return np.array(datos)

def compare_columns(col1, col2, tolerancia=1e-8):
    """Compare xyz columns between two files"""
    if col1.shape != col2.shape:
        print("The files have different sizes.")
        return False
    if np.allclose(col1[:, :3], col2[:, :3], atol=tolerancia):
        print("The columns are the same between the files.")
        return True
    else:
        diferencias = np.where(~np.isclose(col1[:, :3], col2[:, :3], atol=tolerancia))
        print("Differences found in lines:", diferencias[0] + 1)
        return False
    
def filter_values(datos1, datos2, ztol):
    """Filtering data with density close to 1e-3."""
    # Filter value
    filtro = (0.0008 <= datos1[:, 3]) & (datos1[:, 3] <= 0.0012) & (datos1[:, 2] >= ztol)
    #filtro = (0.0001 <= datos1[:, 3]) & (datos1[:, 3] <= 0.0003) & (datos1[:, 2] >= ztol)
    datos_filtrados = datos1[filtro]

    # choosing data closer to 1e-3
    xy_coords = {}
    for x, y, z, valor in datos_filtrados:
        clave = (x, y)
        #if clave not in xy_coords or abs(valor - 2e-4) < abs(xy_coords[clave][1] - 2e-4):
        if clave not in xy_coords or abs(valor - 1e-3) < abs(xy_coords[clave][1] - 1e-3):
            xy_coords[clave] = (z, valor)

    # choosing values in MODELPOT file
    coordenadas_valores_archivo2 = {(x, y, z): valor for x, y, z, valor in datos2}
    resultado = [(x, y, z, coordenadas_valores_archivo2.get((x, y, z), None))
                 for (x, y), (z, valor) in xy_coords.items() if (x, y, z) in coordenadas_valores_archivo2]

    print(f"Number of data filtered in MODELPOT file: {len(resultado)}")
    return resultado

def filter_values_by_z(datos1, datos2, z_target, z_tolerance):
    """Filtering data with z-coordinate close to a target value."""
    # Filter value
    filtro = (z_target - z_tolerance <= datos1[:, 2]) & (datos1[:, 2] <= z_target + z_tolerance)
    datos_filtrados = datos1[filtro]

    # choosing data closer to z_target
    xy_coords = {}
    for x, y, z, valor in datos_filtrados:
        clave = (x, y)
        if clave not in xy_coords or abs(z - z_target) < abs(xy_coords[clave][0] - z_target):
            xy_coords[clave] = (z, valor)

    # choosing values in MODELPOT file
    coordenadas_valores_archivo2 = {(x, y, z): valor for x, y, z, valor in datos2}
    resultado = [(x, y, z, coordenadas_valores_archivo2.get((x, y, z), None))
                 for (x, y), (z, valor) in xy_coords.items() if (x, y, z) in coordenadas_valores_archivo2]

    print(f"Number of data filtered at z ≈ {z_target}: {len(resultado)}")
    return resultado

###########################################
def read_high_symmetry_points(filename):
    """Lee los puntos de alta simetría desde un archivo."""
    points = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():  # Ignorar líneas vacías
                parts = line.split()
                name = parts[0]
                coords = list(map(float, parts[1:]))
                points[name] = np.array(coords)
    return points

def generate_k_path(points, path, min_points=30):
    """Genera el camino de puntos k según el camino de alta simetría dado."""
    subpaths = path.split('|')  # Dividir en subcaminos separados por '|'
    k_points = []

    for subpath in subpaths:
        path_labels = subpath.split('-')  # Dividir los puntos en la secuencia

        # Agregar los puntos de alta simetría en la secuencia, sin interpolación
        if len(path_labels) == 1:
            label = path_labels[0]
            k_points.append((*points[label], 1, f"!{label}"))
            continue

        # Calcular distancias de los segmentos en este subcamino
        segment_lengths = {}
        for i in range(len(path_labels) - 1):
            start_label, end_label = path_labels[i], path_labels[i + 1]
            start_point, end_point = points[start_label], points[end_label]
            segment_length = np.linalg.norm(end_point - start_point)
            segment_lengths[(start_label, end_label)] = segment_length

        # Encontrar la distancia más corta
        min_distance = min(segment_lengths.values())

        # Calcular el número total de puntos en función de la distancia más corta
        total_points = {seg: max(min_points, int((length / min_distance) * min_points))
                        for seg, length in segment_lengths.items()}

        # Generar los puntos interpolados
        for i in range(len(path_labels) - 1):
            start_label, end_label = path_labels[i], path_labels[i + 1]
            start_point, end_point = points[start_label], points[end_label]
            num_points = total_points[(start_label, end_label)]

            # Interpolar puntos
            for step in range(num_points):
                t = step / (num_points - 1)  # Interpolación lineal
                interpolated_point = (1 - t) * start_point + t * end_point

                # Si es el primer punto del segmento, agregarlo siempre
                if step == 0:
                    k_points.append((*interpolated_point, 1, f"!{start_label}"))
                # Si es el último punto, agregarlo solo si el siguiente segmento es discontinuo (|) o es el final
                elif step == num_points - 1:
                    # Verificar si el siguiente segmento es discontinuo
                    is_last_segment = (i == len(path_labels) - 2)
                    next_discontinuous = (subpath in subpaths[:-1])  # Si hay más subpaths después, hay un '|'
                    if is_last_segment or next_discontinuous:
                        k_points.append((*interpolated_point, 1, f"!{end_label}"))
                else:
                    k_points.append((*interpolated_point, 1, ""))

    return k_points

def save_k_path(k_points, filename):
    """Guarda el camino de puntos k en un archivo."""
    with open(filename, 'w') as file:
        for point in k_points:
            x, y, z, w, comment = point
            line = f"{x:.8f} {y:.8f} {z:.8f} {w}"
            if comment:
                line += f" {comment}"
            file.write(line + "\n")

def read_band_data(bands_file):
    """Lee los datos de bandas electrónicas, identificando bandas separadas por líneas vacías."""
    bands = []
    current_band = []

    # Leer el archivo línea por línea
    with open(bands_file, 'r') as file:
        for line in file:
            if line.strip():  # Si la línea no está vacía, agregar los valores
                values = list(map(float, line.split()))
                current_band.append(values)
            elif current_band:  # Si la línea está vacía y ya hay datos, iniciar una nueva banda
                bands.append(np.array(current_band))
                current_band = []

    # Agregar la última banda si hay datos
    if current_band:
        bands.append(np.array(current_band))

    return bands

###########################################
#####funciones principales

def Fukui_interpolation(FILE1,FILE2,FILE3,FILE4, dn=None):
    
    #cuando llamemos la funcion hat que antes converit en dn en un array dn=np.array([n1, n2, n3, n4]

    if dn is None:
        dn = np.array([-0.15, -0.1, -0.05, 0.0])
    print(r"\delta N is: ", dn)
    print ("Reading info from files")

    NATOMS, GRID, POINTS, alattvec, blattvec, clattvec, skipelines, maxrows = read_file_info(FILE1)
    
    print("Colecting info from different files")
    CHG1 = process_file(FILE1, skipelines, maxrows, POINTS)
    CHG2 = process_file(FILE2, skipelines, maxrows, POINTS)
    CHG3 = process_file(FILE3, skipelines, maxrows, POINTS)
    CHG4 = process_file(FILE4, skipelines, maxrows, POINTS)
    
    print('CHG1',CHG1)
    print('CHG2',CHG2)
    print('CHG3',CHG3)
    print('CHG4',CHG4)

    print('Making interpolation')

    FUKUI, RSQUARED = calculate_fukui_and_r2(dn, CHG1, CHG2, CHG3, CHG4)

    print('Fukui size', FUKUI.size)
    print('Fukui', FUKUI) 

    final_file= write_fukui_file(FILE1, NATOMS, FUKUI, "CHGCAR_FUKUI.vasp")
    
    print('FUKUI file was written')
    print("")
    return final_file

def fukui_electrodes(FILE0,FILE1,Epsilon):
    
    NATOMS, GRID, POINTS, alattvec, blattvec, clattvec, skipelines, maxrows = read_file_info(FILE1)
    
    NGX,NGY,NGZ = GRID
    
    print("This will take a few seconds.")
    print("")

    LATTCURB, VOL, omega= compute_lattice_parameters(alattvec, blattvec, clattvec)
    GSQU=compute_gsquared(GRID, LATTCURB, 0)        

    df0 = read_data(FILE0,skipelines, maxrows, POINTS)
    df = read_data(FILE1,skipelines, maxrows, POINTS)
                    
    CHG = missing(df,POINTS)
    CHG00 = missing(df0,POINTS)

    CHGtem = reshape_xyz(NGX, NGY, NGZ, CHG)
    CHGtem00 = reshape_xyz(NGX, NGY, NGZ, CHG00)

    CHG = CHGtem/omega
    CHG00 = CHGtem00/omega

    PLANAR_CHG00 = compute_planar_average_nz(CHG00, NGX, NGY, NGZ, 'z')

    CHGG = np.fft.fftn(CHG, norm='ortho')    
    LOCPOTG = 4*np.pi*np.multiply(CHGG,GSQU)
    LOCPOT = np.fft.ifftn(LOCPOTG,norm='ortho').real * 27.2114
    
    # Set number
    cota = 0.0001
    cota_ls = 0.02
    
    valor_cercano, pos0 = closest_value_position(PLANAR_CHG00, cota)
    val_ls, pos_ls = closest_value_position(PLANAR_CHG00, cota_ls)
    
    LZ = clattvec[2] 
    LS = ((LZ/NGZ)*pos_ls - LZ/2)*2
    
    salida = [item for sublist in func_correction(VOL, Epsilon, LS, LZ, NGZ) for item in np.array(sublist).flatten() ]
    salida = np.array(salida)
    LOCPOT_cor = LOCPOT + salida[np.newaxis, np.newaxis, :]
        

    PLANAR1 = compute_planar_average_nz(LOCPOT, NGX, NGY, NGZ, 'z')
    PLANAR2 = compute_planar_average_nz(LOCPOT_cor, NGX, NGY, NGZ, 'z')
        
    cota2 = PLANAR2[pos0]
    PLANAR2f_2 = [elemento - cota2 if 0 <= (elemento - cota2) else 0 for elemento in PLANAR2] 
    correccion_final = [a - b for a, b in zip(PLANAR2f_2, PLANAR2)]
    correccion_final = np.array(correccion_final)
    LOCPOT_cor2 = LOCPOT_cor + correccion_final[np.newaxis, np.newaxis, :]

    PLANAR3= compute_planar_average_nz(LOCPOT_cor2, NGX, NGY, NGZ, 'z')
      
    PLANAR1 = [-x for x in PLANAR1]
    PLANAR3 = [-x for x in PLANAR3]

    z_axis = np.linspace(0, LZ, NGZ)
    plot_planar_average(z_axis, PLANAR1, PLANAR3)

    LOCPOTtem = convert_locpot(LOCPOT_cor2, NGX, NGY, NGZ)

    final_file = write_fukui_file(FILE1, NATOMS, LOCPOTtem, 'FUKUI.LOCPOT')

    return final_file
  
def fukui_SCPC(FILE0,FILE1,FILE2,c):
    
    NATOMS, GRID, POINTS, alattvec, blattvec, clattvec, skipelines, maxrows = read_file_info(FILE1)
    
    NGX,NGY,NGZ = GRID
    
    print("This will take a few seconds.")
    print("")

    LATTCURB, VOL, omega= compute_lattice_parameters(alattvec, blattvec, clattvec)
    GSQU = compute_gsquared(GRID, LATTCURB, 0)
    
    df0 = read_data(FILE0, skipelines, maxrows, POINTS)
    df = read_data(FILE1, skipelines, maxrows, POINTS)
    
    CHG = missing(df,POINTS)
    CHG00 = missing(df0,POINTS)
    CHGtem = reshape_xyz(NGX, NGY, NGZ, CHG)
    CHGtem00 = reshape_xyz(NGX, NGY, NGZ, CHG00)
    CHG = CHGtem/omega
    CHG00 = CHGtem00/omega

    PLANAR_CHG00 = compute_planar_average_nz(CHG00, NGX, NGY, NGZ, 'z')
    
    cota = 0.0001
    
    valor_cercano, pos0 = closest_value_position(PLANAR_CHG00, cota)
    
    CHGG = np.fft.fftn(CHG, norm='ortho')    
    LOCPOTG = 4*np.pi*np.multiply(CHGG,GSQU)
    LOCPOT = np.fft.ifftn(LOCPOTG,norm='ortho')    
    LOCPOT = 27.2114*LOCPOT.real 
    
    LZ = clattvec[2] 
    z_axis = np.linspace(0, LZ, NGZ)

    PLANAR1 = compute_planar_average_nz(LOCPOT, NGX, NGY, NGZ, 'z')
    
    # value1 = np.array(PLANAR1)
    
    z2, value2 = extract_data_from_file(FILE2, marker='#z-vcor.dat')
    interp2 = interp1d(z2, value2, kind='cubic', fill_value="extrapolate")
    value2_interp = interp2(z_axis)
    
    print('value2_interp=',value2_interp)

    LOCPOT_cor = LOCPOT.copy()
    for i in range(NGZ):
        LOCPOT_cor[:,:,i] = LOCPOT_cor[:,:,i] - value2_interp[i]*(c)

    PLANAR2 = compute_planar_average_nz(LOCPOT_cor, NGX, NGY, NGZ, 'z')
    cota2 = PLANAR2[pos0]
    PLANAR2f_2 = [elemento - cota2 if 0 <= (elemento - cota2) else 0 for elemento in PLANAR2]
    correccion_final = [a - b for a, b in zip(PLANAR2f_2, PLANAR2)]

    
    LOCPOT_cor2 = LOCPOT_cor.copy()
    for i in range(NGZ):
                LOCPOT_cor2[:,:,i] = LOCPOT_cor2[:,:,i] + correccion_final[i]

    PLANAR3= compute_planar_average_nz(LOCPOT_cor2, NGX, NGY, NGZ, 'z')

    PLANAR1 = [-x for x in PLANAR1]
    PLANAR3 = [-x for x in PLANAR3]
    value2_interp = [-x for x in value2_interp]

    z_axis = np.linspace(0, LZ, NGZ)

    plot_planar_average(z_axis, PLANAR1, PLANAR3)

    LOCPOTtem = convert_locpot(LOCPOT_cor2, NGX, NGY, NGZ)

    final_file = write_fukui_file(FILE1, NATOMS, LOCPOTtem, 'FUKUI.LOCPOT')
    return final_file
    
def lineal_operation(FILE1,FILE2,c1,c2,c3):
    
    print ("FILE1: ",FILE1)
    NATOMS, GRID, POINTS, alattvec, blattvec, clattvec, skipelines, maxrows = read_file_info(FILE1)
    
    print ("FILE2: ",FILE2)
    NATOMS2, GRID2, POINTS2, alattvec2, blattvec2, clattvec2, SKIPELINES2, MAXROWS2 = read_file_info(FILE2)
    
    if POINTS != POINTS2:
        print("The files have different systems")
        return
     
    CHG1 = process_file(FILE1, skipelines, maxrows, POINTS)
    print('CHG1',CHG1)

    CHG2 = process_file(FILE2, SKIPELINES2, MAXROWS2, POINTS2)
    print('CHG2',CHG2)
    

    CHGSUM = c1 * CHG1 + c2 * CHG2 + c3
    
    final_file = write_fukui_file(FILE1, NATOMS, CHGSUM, 'CHGCARSUM')
    return final_file
    
def planar_average(FILE1, type_file, axis='z'):
    
    print ("FILE1: ",FILE1)
    NATOMS, GRID, POINTS, alattvec, blattvec, clattvec, skipelines, maxrows = read_file_info(FILE1)
    NGX,NGY,NGZ = GRID
    print(alattvec,blattvec,clattvec)
    LATTCURB, VOL, omega= compute_lattice_parameters(alattvec, blattvec, clattvec)
    LZ = clattvec[2]
    LY = blattvec[1]
    LX = alattvec[0]

    df = read_data(FILE1, skipelines, maxrows, POINTS)
    CHG = df.to_numpy().flatten('C')[0:POINTS]
    GHGtem =  reshape_xyz(NGX, NGY, NGZ, CHG)


    if type_file == 'CHGCAR':
        CHG = GHGtem/omega
        TITLE_plot = 'Planar Average Electron Density'
        TITLE_dat = 'PA_ED.dat'
        y_label = r'$\rho(r) \ (a_{0}^{-3})$'
        name_fig = 'PA_ED.svg'
    
    if type_file == 'LOCPOT':
        CHG = GHGtem
        TITLE_plot = 'Planar Average Electrostatic Potential'
        TITLE_dat = 'PA_EP.dat'
        y_label= r'V(r) (eV)'
        name_fig = 'PA_EP.svg'

 

    if axis == 'z':
        xlabel = r'Z-direction ($\AA$)'
        n_axis = np.linspace(0, LZ, NGZ)
    if axis == 'y':
        xlabel = r'Y-direction ($\AA$)'
        n_axis = np.linspace(0, LY, NGY)
    if axis == 'x':
        xlabel = r'x-direction ($\AA$)'
        n_axis = np.linspace(0, LX, NGX)
    
    print('Preparing ',axis, ' axis data')
    print(NGX,NGY,NGZ)
    print(n_axis)
    print(np.linspace(0, LX, NGX))

    AVEGPOTZcorr = compute_planar_average_nz(CHG, NGX, NGY, NGZ, axis)
    data_to_save = np.column_stack((n_axis, AVEGPOTZcorr))
    

    fmt = '%20.6f %20.5e'
    header = f'{"n_axis".rjust(20)}{"Planar_Avg".rjust(20)}'

    np.savetxt(TITLE_dat, data_to_save, fmt=fmt, header=header, comments='')    
    plt.clf()
    plt.plot(n_axis,AVEGPOTZcorr,linewidth=2)
    plt.title(TITLE_plot,fontsize=18)
    plt.ylabel(y_label,fontsize=12.5)
    plt.xlabel(xlabel,fontsize=12.5)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.tight_layout()
    plt.show()
    plt.savefig(name_fig)

def XYZ_value(FILE1,type_file, c1=1,c2=0,plot=False):
    output_file = f"XYZ_{type_file}.dat"
    
    print ("File to convert to xyz format: ",FILE1)
    NATOMS, GRID, POINTS, alattvec, blattvec, clattvec, skipelines, maxrows = read_file_info(FILE1)
    NGX,NGY,NGZ = GRID
    
    LATTCURB, VOL, omega = compute_lattice_parameters(alattvec, blattvec, clattvec)

    LATTCURADO = (1.0) * np.dstack([alattvec, blattvec, clattvec])
    LATTCURADO = LATTCURADO[0]
    
    LZ = clattvec[2]

    df = read_data(FILE1, skipelines, maxrows, POINTS)
    CHG1 = df.to_numpy().flatten('C')[0:POINTS]     

    #CHG1 = process_file(FILE1, skipelines, maxrows, POINTS)
    GHGtem = reshape_xyz(NGX,NGY,NGZ, CHG1)

    if type_file == 'CHGCAR':
        CHG = GHGtem/omega
    if type_file == 'LOCPOT': 
        CHG = GHGtem*(c1)
    
        AVEGPOTZcorr = compute_planar_average_nz(CHG, NGX, NGY, NGZ, 'z')
        first_value_AVEGPOTZcorr = AVEGPOTZcorr[0]
     
        if c2 == 1:
            CHG = CHG - first_value_AVEGPOTZcorr

        AVEGPOTZ = compute_planar_average_nz(CHG, NGX, NGY, NGZ, 'z')

        if plot:
            plt.plot(AVEGPOTZ)
            plt.xlabel('Z-direction')
            plt.ylabel('Planar Average of LOCPOT')
            plt.title('Electrostatic Potential')
            plt.savefig('PA_MEP.svg', format='svg', bbox_inches='tight')
            plt.show()

    
    final_file, rho = write_xyz(CHG, NGX, NGY, NGZ, LATTCURADO, omega, output_file)
    return final_file

def Perturbative_point(FILE1,FILE2,q,N):
    
    print ("FILE1: ",FILE1)
    NATOMS, GRID, POINTS, alattvec, blattvec, clattvec, skipelines, maxrows = read_file_info(FILE1)
    NGX,NGY,NGZ = GRID

    print ("FILE2: ",FILE2)
    NATOMS2, GRID2, POINTS2, alattvec2, blattvec2, clattvec2, SKIPELINES2, MAXROWS2 = read_file_info(FILE2)
    
    if POINTS != POINTS2:
        print("The files have different grids")
        return
    
    
    df = read_data(FILE1, skipelines, maxrows, POINTS)
    CHG1 = df.to_numpy().flatten('C')[0:POINTS]  
    #CHG1 = process_file(FILE1,skipelines,maxrows, POINTS)
    CHG1 = CHG1.astype(np.float64)
    print('CHG1',CHG1)
    
    df = read_data(FILE2, skipelines, maxrows, POINTS)
    CHG2 = df.to_numpy().flatten('C')[0:POINTS]  
    #CHG2 = process_file(FILE2,SKIPELINES2, MAXROWS2, POINTS)
    CHG2 = CHG2.astype(np.float64)
    print('CHG2',CHG2)

    
    print("")
    print("Just a few seconds.")

    CHGtem1 = reshape_xyz(NGX, NGY, NGZ, CHG1)
    CHGtem2 = reshape_xyz(NGX, NGY, NGZ, CHG2)


    CHGSUM = q*CHGtem1*(-1) - N*q*CHGtem2*(-1)
    

    AVEGPOTZcorr = compute_planar_average_nz(CHGSUM, NGX, NGY, NGZ, 'z')
    ##correción y alineación al vacío
    first_value_AVEGPOTZcorr = AVEGPOTZcorr[0]
    
    
    CHG = CHGSUM - first_value_AVEGPOTZcorr

    LOCPOTtem =  np.zeros(NGX*NGY*NGZ)
    for N1 in range(NGZ):
        for N2 in range(NGY):
            for N3 in range(NGX):
                LOCPOTtem[N1*NGX*NGY + N2*NGX + N3] = CHG[N3,N2,N1]


    final_file = write_fukui_file(FILE2, NATOMS, LOCPOTtem, 'MODELPOT.LOCPOT') 
    return final_file

def min_max(FILE1, extrema_point ):
    NATOMS, GRID, POINTS, alattvec, blattvec, clattvec, skipelines, maxrows = read_file_info(FILE1)
    POT = process_file(FILE1, skipelines, maxrows, POINTS)
    POT = np.array(POT).astype(float)
    POT = POT.reshape(GRID[2], GRID[1], GRID[0])
    coords, values = detect_local_extrema_3D(POT, order=2, extrema_type=extrema_point)
    m_points = calcular_xyz_val(coords, values, GRID, alattvec, blattvec, clattvec)
    name_file = f"{extrema_point}.txt"
    write_formatted_data_localm(name_file, m_points)

def visual_modelpot(file1, file2):
#def visual_modelpot(file1, file2, z_target=27.70, z_tolerance=0.1):
    # Files to compare
    archivo1 = XYZ_value(file1, 'CHGCAR')
    archivo2 = XYZ_value(file2, 'LOCPOT')

    print('Reading: ', archivo1,'and', archivo2, 'to plot')
    NATOMS, GRID, POINTS, alattvec, blattvec, clattvec, skipelines, maxrows = read_file_info(file1)
    NGX,NGY,NGZ = GRID
    ztol = np.linalg.norm(clattvec)*0.5

    datos1 = read_xyzval(archivo1)
    datos2 = read_xyzval(archivo2)

    compare_columns(datos1, datos2)

    filtrados = filter_values(datos1, datos2, ztol)
    #filtrados = filter_values_by_z(datos1, datos2, z_target, z_tolerance)

    # Guardar los datos filtrados en un archivo
    with open("datos_filtrados.txt", "w") as f:
        f.write("# x, y, z, valor\n")
        for (x, y, z, valor) in filtrados:
            f.write(f"{x:>12.4f} {y:>12.4f} {z:>12.4f} {valor:>20.8e}\n")
    
    # Values for heatmap
    x_vals = [item[0] for item in filtrados]
    y_vals = [item[1] for item in filtrados]
    valores = [item[3] for item in filtrados]

    heatmap_data, x_edges, y_edges = np.histogram2d(x_vals, y_vals, bins=(NGX, NGY), weights=valores)

    x_range = x_edges[-1] - x_edges[0]
    y_range = y_edges[-1] - y_edges[0]

    plt.figure(figsize=(8, 6))
    plt.imshow(heatmap_data.T, cmap='jet_r', origin='lower',
               extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
               aspect=x_range/y_range)  # Escalar ejes según el rango


    cbar = plt.colorbar()
    cbar.set_label(r"$\Delta U_{int}$", fontsize=14, fontweight='bold') 
    cbar.ax.tick_params(labelsize=12) 


    #plt.title(r"$\Delta U_{int}$ at  $\rho =10^{-4}$ $a_{0}^{-3}$", size=16, family='sans-serif')
    plt.title(r"$\Delta U_{int}$ at  $\rho =10^{-3}$ $a_{0}^{-3}$", size=16, family='sans-serif')
    #plt.title(r"$\Delta U_{int}$ at  $Z = 28.0$", size=16, family='sans-serif')
    plt.xlabel(r"X (Å)",size=14, family='sans-serif')
    plt.ylabel(r"Y (Å)", size=14,family='sans-serif')

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.savefig("heatmap_MODELPOT.png")
    plt.show()

###########################################
def electronic_band_structure(FILE1, FILE2, fermi_energy, y_limits=None):
    """Grafica la estructura de bandas electrónicas con opción de definir límites en el eje Y sin afectar el eje X."""

    # Leer los datos de bandas
    bands = read_band_data(FILE1)

    # Convertir energías de Ry a eV y ajustar por energía de Fermi
    for i in range(len(bands)):
        bands[i][:, 1] = (bands[i][:, 1]) * 13.6057 - fermi_energy

    # Tomar el eje x de la primera banda como referencia para las etiquetas de alta simetría
    reference_x = bands[0][:, 0]

    # Cargar los datos del archivo k-path.dat
    kpath_data = np.loadtxt(FILE2, dtype=str)

    # Extraer índices y etiquetas de alta simetría
    k_indices = kpath_data[:, 0].astype(int) - 1  # Convertir a índice 0-based
    k_symbols = np.where(kpath_data[:, 1] == "G", r"$\Gamma$", kpath_data[:, 1])

    # Obtener las posiciones de los puntos de alta simetría en el eje x
    new_k_positions = reference_x[k_indices]

    # Crear la figura
    fig, ax = plt.subplots(figsize=(8, 6))

    # Graficar cada banda de forma independiente
    for band in bands:
        ax.plot(band[:, 0], band[:, 1], color="steelblue", lw=2.5)

    # Añadir líneas verticales en los puntos de alta simetría
    for k_pos in new_k_positions:
        ax.axvline(x=k_pos, color="black", linestyle="-", linewidth=1.0, zorder=0)

    # **⚠️ Guardar los límites originales del eje X**
    original_x_limits = (min(reference_x), max(reference_x))  # Tomar los valores extremos de la primera banda

    # Configuración de los ejes
    ax.set_ylabel("Energía (eV)", fontsize=14)
    ax.set_xticks(new_k_positions)
    ax.set_xticklabels(k_symbols, fontsize=12)  # Etiquetas de alta simetría alineadas correctamente
    ax.tick_params(axis='y', labelsize=12)

    # Agregar una línea en la energía de Fermi
    ax.axhline(y=0, color='red', linestyle='--', linewidth=1.5)

    # **⚠️ Aplicar los límites del eje Y si están definidos**
    if y_limits is not None and len(y_limits) == 2:
        ax.set_ylim(y_limits[0], y_limits[1])  # Forzar límites en Y
        print(f"Applying Y-axis limits: {y_limits[0]} to {y_limits[1]}")

    # **⚠️ Forzar los límites del eje X nuevamente**
    ax.set_xlim(original_x_limits)  # Restaurar los límites originales del eje X

    # **Desactivar ajuste automático de Matplotlib**
    ax.set_autoscale_on(False)

    # **🔹 Desactivar márgenes innecesarios**
    ax.margins(x=0)  

    # **🔹 Evitar que Matplotlib ajuste el gráfico**
    fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

    # Guardar el gráfico como archivo SVG
    output_file = "Bands.svg"
    plt.savefig(output_file, format="svg")

    print(f"The plot has been saved as {output_file}")

def phonon_dispersion_plot(file1, file2):
    # Cargar los datos de frecuencias
    freq_file = file1
    data = np.loadtxt(freq_file)

    # Separar las columnas
    x = np.arange(len(data))  # Índices de las filas (0 a n-1)
    y_columns = data[:, 1:]  # Resto de las columnas para graficar

    # Factor de escala para frecuencias
    scale_factor = 0.124
    y_columns *= scale_factor

    # Cargar los datos del archivo k-path.dat
    kpath_file = file2
    kpath_data = np.loadtxt(kpath_file, dtype=str)

    # Convertir etiquetas de posición y reemplazar "G" por "\Gamma"
    k_positions = kpath_data[:, 0].astype(int) - 1  # Convertir a índices (restar 1 para índice 0-based)
    k_symbols = np.where(kpath_data[:, 1] == "G", r"$\Gamma$", kpath_data[:, 1])  # Reemplazar "G" por "\Gamma"

    # Combinar etiquetas consecutivas
    new_k_positions = []
    new_k_symbols = []
    i = 0
    while i < len(k_positions):
        if i < len(k_positions) - 1 and k_positions[i] + 1 == k_positions[i + 1]:
            combined_label = f"{k_symbols[i]}/{k_symbols[i+1]}"
            new_k_positions.append((k_positions[i] + k_positions[i + 1]) / 2)
            new_k_symbols.append(combined_label)
            i += 2  # Saltar el siguiente punto porque ya fue combinado
        else:
            new_k_positions.append(k_positions[i])
            new_k_symbols.append(k_symbols[i])
            i += 1
    # Crear el gráfico
    plt.figure(figsize=(8, 6))

    # Graficar todas las curvas con color azul claro y líneas más gruesas
    for col in range(y_columns.shape[1]):
        plt.plot(x, y_columns[:, col], color="steelblue", lw=2.5)

    # Añadir líneas verticales en los puntos de alta simetría
    for k_pos in new_k_positions:
        plt.axvline(x=k_pos, color="black", linestyle="-", linewidth=1.0, zorder=0)
   
    # Configuración de los ejes
    plt.ylabel(r"$\omega$ (meV)", fontsize=14)  # Tamaño de etiqueta del eje y
    plt.xticks(new_k_positions, new_k_symbols, fontsize=12)  # Etiquetas personalizadas con mayor tamaño
    plt.yticks(fontsize=12)  # Tamaño de números en el eje y
  
    # Eliminar las pequeñas líneas (ticks menores) del eje x
    plt.tick_params(axis='x', which='both', bottom=False)  # Sin ticks menores ni mayores en el eje x

    # Ajustar límites del eje y
    max_y = np.max(y_columns)  # Valor máximo de las frecuencias
    plt.ylim(0, max_y + 10)  # Eje y desde 0 con margen superior de 10

    # Ajustar márgenes del eje x
    plt.margins(x=0)  # Sin margen blanco en el eje x

    # Opciones adicionales
    plt.tight_layout()

    # Guardar el gráfico como archivo SVG
    output_file = "Phonon_dipersion_plot.svg"
    plt.savefig(output_file, format="svg")

    # Mostrar mensaje de confirmación
    print(f"The plot has been saved as {output_file}")


def compare_phonon_dispersion_two(file1, file2, kpath_file, label1="Calculation 1", label2="Calculation 2", title="Phonon Dispersion Comparison"):
    """
    Compara dos cálculos de dispersión de fonones en un solo gráfico.

    Parámetros:
    - file1: Archivo con los datos de frecuencias del primer cálculo.
    - file2: Archivo con los datos de frecuencias del segundo cálculo.
    - kpath_file: Archivo con las posiciones de alta simetría.
    - label1: Etiqueta personalizada para el primer cálculo.
    - label2: Etiqueta personalizada para el segundo cálculo.
    """

    # Cargar los datos de frecuencias
    data1 = np.loadtxt(file1)
    data2 = np.loadtxt(file2)

    # Separar las columnas (índices en X y frecuencias en Y)
    x1 = np.arange(len(data1))
    y1_columns = data1[:, 1:]  # Ignorar la primera columna (puntos)
    
    x2 = np.arange(len(data2))
    y2_columns = data2[:, 1:]

    # Factor de escala para frecuencias (meV)
    scale_factor = 0.124
    y1_columns *= scale_factor
    y2_columns *= scale_factor

    # Cargar los datos del archivo k-path.dat
    kpath_data = np.loadtxt(kpath_file, dtype=str)

    # Convertir etiquetas de posición y reemplazar "G" por "\Gamma"
    k_positions = kpath_data[:, 0].astype(int) - 1  # Convertir a índices (restar 1 para índice 0-based)
    k_symbols = np.where(kpath_data[:, 1] == "G", r"$\Gamma$", kpath_data[:, 1])

    # Crear el gráfico
    plt.figure(figsize=(8, 6))

    # Graficar el primer cálculo (azul sólido)
    for col in range(y1_columns.shape[1]):
        plt.plot(x1, y1_columns[:, col], color="blue", lw=2, alpha=0.6, label=label1 if col == 0 else "")

    # Graficar el segundo cálculo (rojo punteado)
    for col in range(y2_columns.shape[1]):
        plt.plot(x2, y2_columns[:, col], color="red", lw=2, alpha=0.6, linestyle="dashed", label=label2 if col == 0 else "")

    # Añadir líneas verticales en los puntos de alta simetría
    for k_pos in k_positions:
        plt.axvline(x=k_pos, color="black", linestyle="-", linewidth=1.0, zorder=0)
    
    # 📌 Agregar título personalizado
    plt.title(title, fontsize=16, pad=20)

    # Configuración de los ejes
    plt.ylabel(r"$\omega$ (meV)", fontsize=14)
    plt.xticks(k_positions, k_symbols, fontsize=12)
    plt.yticks(fontsize=12)

    # Eliminar ticks menores en el eje X
    plt.tick_params(axis='x', which='both', bottom=False)

    # Ajustar límites del eje Y
    max_y = max(np.max(y1_columns), np.max(y2_columns))  # Máximo entre los dos cálculos
    plt.ylim(0, max_y + 10)

    # Ajustar márgenes del eje X
    plt.margins(x=0)

    # Añadir leyenda con etiquetas personalizadas
    #plt.legend(fontsize=12, loc="upper right")
    plt.legend(fontsize=12, loc="lower right")
    #plt.legend(fontsize=12, loc="upper left", bbox_to_anchor=(1, 1))

    # Ajustar el diseño para que la leyenda no corte el gráfico
    #plt.tight_layout(rect=[0, 0, 0.85, 1])  # Ajusta el espacio para la leyenda fuera del área del gráfico    

    # Opciones adicionales
    plt.tight_layout()

    # Guardar el gráfico como archivo SVG
    output_file = "Phonon_comparison_two.svg"
    plt.savefig(output_file, format="svg")

    # Mostrar mensaje de confirmación
    print(f"The comparison plot has been saved as {output_file}")


def compare_phonon_dispersion_three(FILE1, FILE2, FILE3, FILE4, label1="Calculation 1", label2="Calculation 2", label3="Calculation 3", title="Phonon Dispersion Comparison"):
    """
    Compara tres cálculos de dispersión de fonones en un solo gráfico.

    Parámetros:
    - file1, file2, file3: Archivos con los datos de frecuencias.
    - file4: Archivo con las posiciones de alta simetría.
    - label1, label2, label3: Etiquetas para cada cálculo.
    - title: Título del gráfico.
    """

    # Cargar los datos de frecuencias
    data1 = np.loadtxt(FILE1)
    data2 = np.loadtxt(FILE2)
    data3 = np.loadtxt(FILE3)

    # Separar las columnas (índices en X y frecuencias en Y)
    x1 = np.arange(len(data1))
    y1_columns = data1[:, 1:]
    
    x2 = np.arange(len(data2))
    y2_columns = data2[:, 1:]

    x3 = np.arange(len(data3))
    y3_columns = data3[:, 1:]

    # Factor de escala para frecuencias (meV)
    scale_factor = 0.124
    y1_columns *= scale_factor
    y2_columns *= scale_factor
    y3_columns *= scale_factor

    # Cargar los datos del archivo k-path.dat
    kpath_data = np.loadtxt(FILE4, dtype=str)

    # Convertir etiquetas de posición y reemplazar "G" por "\Gamma"
    k_positions = kpath_data[:, 0].astype(int) - 1
    k_symbols = np.where(kpath_data[:, 1] == "G", r"$\Gamma$", kpath_data[:, 1])

    # Crear el gráfico
    fig, ax = plt.subplots(figsize=(8, 6))

    # Graficar las tres curvas con diferentes colores y estilos
    for col in range(y1_columns.shape[1]):
        ax.plot(x1, y1_columns[:, col], color="blue", lw=2, alpha=0.6, label=label1 if col == 0 else "")

    for col in range(y2_columns.shape[1]):
        ax.plot(x2, y2_columns[:, col], color="red", lw=2, alpha=0.6, linestyle="dashed", label=label2 if col == 0 else "")

    for col in range(y3_columns.shape[1]):
        ax.plot(x3, y3_columns[:, col], color="green", lw=2, alpha=0.6, linestyle="dashdot", label=label3 if col == 0 else "")

    # Añadir líneas verticales en los puntos de alta simetría
    for k_pos in k_positions:
        ax.axvline(x=k_pos, color="black", linestyle="-", linewidth=1.0, zorder=0)

    # 📌 Agregar título personalizado
    plt.title(title, fontsize=16, pad=20)

    # Configuración de los ejes
    ax.set_ylabel(r"$\omega$ (meV)", fontsize=14)
    ax.set_xticks(k_positions)
    ax.set_xticklabels(k_symbols, fontsize=12)
    ax.tick_params(axis='y', labelsize=12)

    # Agregar una línea en la energía de referencia
    ax.axhline(y=0, color='red', linestyle='--', linewidth=1.5)

    # Ajustar límites del eje Y
    max_y = max(np.max(y1_columns), np.max(y2_columns), np.max(y3_columns))
    ax.set_ylim(0, max_y + 10)

    # Ajustar márgenes del eje X
    ax.margins(x=0)

    # Añadir leyenda con etiquetas personalizadas
    #plt.legend(fontsize=12, loc="upper right")
    plt.legend(fontsize=12, loc="lower right")
    #plt.legend(fontsize=12, loc="upper left", bbox_to_anchor=(1, 1))

    # Ajustar el diseño para que la leyenda no corte el gráfico
    #plt.tight_layout(rect=[0, 0, 0.85, 1])  # Ajusta el espacio para la leyenda fuera del área del gráfico    

    # Opciones adicionales
    plt.tight_layout()

    # Guardar el gráfico como archivo SVG
    output_file = "Phonon_comparison_three.svg"
    plt.savefig(output_file, format="svg")

    # Mostrar mensaje de confirmación
    print(f"The comparison plot has been saved as {output_file}")

def k_path(FILE1, path):
    """Genera y guarda el camino de puntos k basado en puntos de alta simetría."""
    # Leer los puntos de alta simetría
    high_symmetry_points = read_high_symmetry_points(FILE1)

    # Generar el camino k
    k_points = generate_k_path(high_symmetry_points, path)

    # Guardar los puntos k en un archivo
    output_filename = "k_path.txt"
    save_k_path(k_points, output_filename)

    print(f"The points have been saved in '{output_filename}'")

def plot_superconductor_iso_gap(file):
    """
    Lee un archivo con temperaturas (K) y gaps superconductores (meV) y genera un gráfico de puntos.

    Parámetros:
    - file: Nombre del archivo con los datos.
    """

    try:
        # Cargar los datos desde el archivo
        data = np.loadtxt(file)

        # Verificar que el archivo no está vacío
        if data.size == 0:
            print("Error: The file is empty or has incorrect formatting.")
            return

    except Exception as e:
        print(f"Error loading file: {e}")
        return  # Salir de la función si hay un error

    # Separar columnas
    temperature = data[:, 0]  # Temperatura en K
    gap = data[:, 1]  # Gap superconductor en meV

    # Crear figura
    fig, ax = plt.subplots(figsize=(6, 4))

    # Graficar solo puntos (cuadrados azules)
    ax.scatter(temperature, gap, color="blue", marker="s", s=50, edgecolors="blue", linewidth=0.6)

    # Configuración de ejes con formato LaTeX
    ax.set_xlabel(r"Temperature (K)", fontsize=14)
    ax.set_ylabel(r"$\Delta_{0}$ (meV)", fontsize=14)

    # Ajustar límites automáticamente
    ax.set_xlim(0, max(temperature) + 20)
    ax.set_ylim(5, max(gap) + 5)

    # Ajustar ticks del eje Y y X
    ax.tick_params(axis='both', labelsize=12)

    # Ajustar diseño
    plt.tight_layout()

    # Guardar gráfico
    output_file = "Superconducting_Gap_vs_Temperature.svg"
    plt.savefig(output_file, format="svg")

    # Mostrar mensaje de confirmación
    print(f"The plot has been saved as {output_file}")

  
############################################################
#                 Main Menu                                #
############################################################

def main_menu():
        while True:

            #******Print Header*****
            header = """

            8
            8
       888  8   88
      8   8 8  8  8
     8    8 8 8    8
           888
            8
            8
            8
            8
           888
               7 
               7 
               777 
"""

            print(header)
            print("CMPkit -- A useful tool for Condensed Matter Physics.")
            print("Version 1.0, release date: March-2025.")
            print("Developers: Nicolas F. Barrera.")
            print("            LCT, Sorbonne Universite.")
            print("")
            print("")
            print("")
            print("               ***** Main Menu ***** ")
            print("1 -- Electronic Band Structure.")
            print("2 -- Phonon Dispersion plots.")
            print("3 -- Interpolation of electron-phonon matrix elements plots.")
            print("4 -- Tc via Eliashberg equations plots.")
            print("5 -- Deformation potential.")
            print("6 -- Exit.")
            print("")
            option = input("Choose an option: ")

            if option == "1":
                print("\nChose option 1: Electronic Band Structures.\n")
                print("11 Generation of k-path from high-symmetry points.")
                print("12 Bands Structure plot.")

                option1 = input("Choose an option: ")
                if option1 == "11":
                    print("\nName of high-symmetry points file.\n")
                    FILE1 = input("Enter name of file: ")
                    path = input("Ingrese el camino de puntos de alta simetría (por ejemplo, T-U|X-L): ").strip()                    
                    print("\nThis will take a few seconds.\n")

                    k_path(FILE1,path)

                    continue
    
                if option1 == "12":
                    print("Name of bands structures file with extension out.gnu")
                    FILE1 = input("Enter name of file: ")
                    print("\nName of high-symmetry points file.\n")
                    FILE2 = input("Enter name of file: ")
                    print("Enter the Fermi energy: ")
                    fermi_energy = float(input("Ef (eV): "))
                    print("\nThis will take a few seconds.\n")
                    # Preguntar si el usuario quiere definir límites en el eje Y
                    print("\nYou can define y-axis limits or press Enter to set them automatically.")

                    # Preguntar si el usuario quiere definir límites en el eje Y
                    print("\nYou can define y-axis limits or press Enter to set them automatically.")

                    # Inicializar variable de límites
                    y_limits = None  

                    # Capturar el límite inferior
                    y_min_input = input("Enter the lower limit of the y-axis (or press Enter to skip): ").strip()
                    y_max_input = input("Enter the upper limit of the y-axis (or press Enter to skip): ").strip()

                    # Verificar si el usuario ingresó ambos valores
                    if y_min_input and y_max_input:
                        try:
                            y_limits = (float(y_min_input), float(y_max_input))
                            print(f"Applying user-defined Y-axis limits: {y_limits[0]} to {y_limits[1]}")
                        except ValueError:
                             print("⚠️ Invalid input. Using automatic limits.")

                    print("\nThis will take a few seconds...\n")
                    electronic_band_structure(FILE1, FILE2, fermi_energy, y_limits)
                    
                    continue

                else:
                    print("Invalid option. Please select an option from the menu.")

            elif option == "2":
                print("You selected option 2: Phonon Dispersion plots.")
                print("")
                print("")
                print("21 -- Phonon Dispersion plot.")
                print("22 -- Phonon Dispersion Comparison plot (two curves).")
                print("23 -- Phonon Dispersion Comparison plot (three curves).")
                print("")

                option2 = input("Choose an option: ")
                if option2 == "21":
                    print("Name of phonon dipersion file.")
                    FILE1 = input("Enter file name: ")
                    print("\nEnter the high-symmetry k-path file.\n")
                    FILE2 = input("Enter file name: ")
                    phonon_dispersion_plot(FILE1, FILE2)
                
                    continue
                
                if option2 == "22":
                    print("Enter the first phonon dispersion file.")
                    FILE1 = input("Enter file name: ")
                    print("Enter label for the first calculation.")
                    LABEL1 = input("Enter label (e.g., '3x3x3'): ").strip()
                    print("Enter the second phonon dispersion file.")
                    FILE2 = input("Enter file name: ")
                    print("Enter label for the second calculation.")
                    LABEL2 = input("Enter label (e.g., '6x6x6'): ").strip()
                    print("\nEnter the high-symmetry k-path file.\n")
                    FILE3 = input("Enter file name: ")
                    print("\nEnter the title for the plot.\n")
                    TITLE = input("Enter title (e.g., 'Phonon Dispersion Comparison'): ").strip()
                    TITLE = fr"${TITLE}$"
                    
                    compare_phonon_dispersion_two(FILE1, FILE2, FILE3, LABEL1, LABEL2, TITLE)
                    
                    continue       
          
                if option2 == "23":
                    print("Enter the first phonon dispersion file.")
                    FILE1 = input("Enter file name: ")
                    print("Enter label for the first calculation.")
                    label1 = input("Enter label (e.g., '3x3x3'): ").strip()
                    print("Enter the second phonon dispersion file.")
                    FILE2 = input("Enter file name: ")
                    print("Enter label for the second calculation.")
                    label2 = input("Enter label (e.g., '6x6x6'): ").strip()
                    print("Enter the third phonon dispersion file.")
                    FILE3 = input("Enter file name: ")
                    print("Enter label for the third calculation.")
                    label3 = input("Enter label (e.g., '3x3x3'): ").strip()
                    print("\nEnter the high-symmetry k-path file.\n")
                    FILE4 = input("Enter file name: ")
                    print("\nEnter the title for the plot.\n")
                    title = input("Enter title (e.g., 'Phonon Dispersion Comparison'): ").strip()
                    title= fr"${title}$"
                    
                    compare_phonon_dispersion_three(FILE1, FILE2, FILE3, FILE4, label1, label2, label3, title)
                    
                    continue      
            
            elif option == "3":
                print("You selected option 3: Fukui Potential via SCPC")
                print("")
                print("")
                print("31 -- Electrophilic Fukui function v_f^-(r).")
                print("32 -- Nucleophilic Fukui function v_f^+(r).")
                print("")

                option3 = input("Choose an option: ")
                if option3 == "31":
                    print("Name CHGCAR file of charge density of the neutral slab.")
                    FILE0 = input("Enter file name: ")
                    print("\nName CHGCAR file of Fukui function.")
                    FILE1 = input("Enter file name: ")
                    print("\nName SCPC correcton file.")
                    FILE2 = input("Enter file name: ")

                    fukui_SCPC(FILE0,FILE1,FILE2,-1)

                    continue
                
                if option3 == "32":
                    print("Name CHGCAR file of charge density of the neutral slab.")
                    FILE0 = input("Enter file name: ")
                    print("\nName CHGCAR file of Fukui function.")
                    FILE1 = input("Enter file name: ")
                    print("\nName SCPC correcton file.")
                    FILE2 = input("Enter file name: ")

                    fukui_SCPC(FILE0,FILE1,FILE2,1)
                    
                    continue

                else:
                    print("Invalid option. Please select an option from the menu.")

            elif option == "4":
                print("You selected option 4: Tc via Eliashberg equations plots")
                print("")
                print("")
                print("41 -- Isotropic superconducting gap.")
                print("42 -- Anisotropic superconducting gap.")
                print("")
                print("")

                option4 = input("Choose an option: ")
                if option4 == "41":
                    print("Name Isotropic superconducting gap file: ")
                    file = input("Enter name of file 1: ")
                    print("\nThis will take a few seconds.\n")
                    
                    plot_superconductor_iso_gap(file)

                    continue

                if option4 == "42":
                    print("Name CHGCAR or LOCPOT files to subtract: ")
                    FILE1 = input("Enter name of file 1: ")
                    FILE2 = input("Enter name of file 2: ")
                    print("\nThis will take a few seconds.\n")

                    lineal_operation(FILE1,FILE2,1,-1,0)
                    
                    # Change the file name
                    old_filename = "CHGCARSUM"
                    new_filename = "CHGCAR_DIFF"
                    os.rename(old_filename, new_filename)

                    continue

            elif option == "6":
                print("You selected option 6: Goodbye!")
                sys.exit()

            else:
                print("Invalid option. Please select an option from the menu.")

if __name__ == "__main__":
        main_menu()
