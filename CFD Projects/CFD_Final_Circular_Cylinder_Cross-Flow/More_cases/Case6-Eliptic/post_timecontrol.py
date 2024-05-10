import numpy as np
import re
import matplotlib.pyplot as plt

class mesher:
    def __init__(self, data):
        self.Lx = data["Lx"]
        self.Ly = data["Ly"]
        # self.Re = data["Re"]
        self.dxy = data["dy"]

    def parameter_grid(self): # generate large grid 
        Field1 = grid_generator( self.Ly,  self.Lx, self.dxy)
        N_rows, N_cols,a = np.shape(Field1)

        # save location y,x in field[j,i,0], 
        for j in range(0, N_rows):
            for i in range(0,N_cols):
                Field1[j,i,0] , Field1[j,i,1] = (j-1/2)*self.dxy, (i-1/2)*self.dxy
        # print(Field1[:,:,1])
        # save fluid(1)/Solid(0) in field[j,i,2]
        Field1[:,:,2] = i_fluid_detector(Field1[:,:,0], Field1[:,:,1], Field1[:,:,2])
        Field1[:,:,3] = 1 - Field1[:,:,2]

        Field1.flags.writeable = False
        return Field1
    
def grid_generator( y_size,  x_size, grid_size):
        Nx , Ny = int(x_size/grid_size), int(y_size/grid_size)
        Field = np.full((Ny+2, Nx+2, 10), 0 ,dtype = 'float')
        return Field


def i_fluid_detector(field_y, field_x, field):
    N_rows, N_cols = np.shape(field)
    for j in range(0, N_rows):
        for i in range(0, N_cols ):
            result = 0
            result = Ellip_function(field_y[j,i], field_x[j,i])
            if result < 0:
                field[j,i] = 0
            else:
                field[j,i] = 1
    return field

def circle_function(y,x):
    R= 0.25
    x_c = data['Lx']*1/4
    y_c = data['Ly']*1/2
    result = (x-x_c)**2 + (y-y_c)**2 - R**2
    return result


def Ellip_function(y,x):
    a = 0.6
    b = 0.2
    cos = np.sqrt(3)/2
    sin = -1/2
    x_c = data['Lx']*1/4
    y_c = data['Ly']*1/2
    x_new = x-x_c
    y_new = y-y_c
    result = ((x_new*cos+y_new*sin)/a)**2 + ((-x_new*sin + y_new*cos)/b)**2 - 1
    return result


def Circle_outer(iF, iS):
    outer = np.zeros_like(iF)
    outer_w = iF[1:-1,1:-1] - iF[1:-1,:-2]; outer_e = iF[1:-1,1:-1] - iF[1:-1,2:]
    outer_s = iF[1:-1,1:-1] - iF[:-2,1:-1]; outer_n = iF[1:-1,1:-1] - iF[2:,1:-1]
    outer[1:-1,1:-1] = outer_w + outer_e + outer_n + outer_s
    for j in range(0,np.size(outer, 0)):
        for i in range(0,np.size(outer,1)):
            if outer[j,i] !=0: outer[j,i] = 1
    outer -= iS

    for j in range(0,np.size(outer, 0)):
        for i in range(0,np.size(outer,1)):
            if outer[j,i] <0: outer[j,i] = 0
    return outer

def Circle_outer_classify(outer):
    outer_e = np.zeros_like(outer);  outer_w = np.zeros_like(outer)
    outer_n = np.zeros_like(outer);  outer_s = np.zeros_like(outer)
    for j in range(0,np.size(outer, 0)):
        for i in range(0,np.size(outer,1)):
            if outer[j,i]==1 and iS[j,i+1]==1:
                outer_w[j,i] =1
            if outer[j,i]==1 and iS[j,i-1]==1:
                outer_e[j,i] =1
            if outer[j,i]==1 and iS[j+1,i]==1:
                outer_s[j,i] =1
            if outer[j,i]==1 and iS[j-1,i]==1:
                outer_n[j,i] =1
          # test ifluid
    # plt.imshow(outer- outer_s, cmap = 'gray')
    # plt.show() 
    return outer_e, outer_w, outer_n, outer_s



def load_2d_arrays_with_numpy(filename):
    # # filename = f'cyl_output_N=32.0, Re=150.0, t=10.0.txt'
    # filename = f'cyl_output_N=32, Re=150, t=50.txt'
    try:
        with open(filename, 'r') as file:
            content = file.read()
    except FileNotFoundError:
        print("文件未找到,请检查文件路径是否正确")
        return None, None, None

    # 辅助函数来提取并构建二维数组
    def extract_2d_array(array_content):
        array_data = []
        # 匹配方括号内的所有内容,包括可能的多行
        matches = re.findall(r'\[(.*?)]', array_content.replace('\n', ''), re.DOTALL)
        for match in matches:
            # 替换掉可能存在的换行符,然后解析数值
            formatted_match = ' '.join(match.split())
            numbers = np.fromstring(formatted_match, sep=' ')
            array_data.append(numbers)
        return np.array(array_data)

    u_data = re.search(r'Array u:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    v_data = re.search(r'Array v:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    U_data = re.search(r'Array U:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    V_data = re.search(r'Array V:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    u_old_data = re.search(r'Array u_old:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    v_old_data = re.search(r'Array v_old:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    U_old_data = re.search(r'Array U_old:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    V_old_data = re.search(r'Array V_old:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    p_data = re.search(r'Array p:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)

    if not (u_data and v_data and p_data):
        print("未能正确提取u, v, 或 p的数据。请检查数组标记和文件格式。")
        return np.array([]), np.array([]), np.array([])

    u = extract_2d_array(u_data.group(1))
    v = extract_2d_array(v_data.group(1))
    p = extract_2d_array(p_data.group(1))
    U = extract_2d_array(U_data.group(1))
    V = extract_2d_array(V_data.group(1))
    u_old = extract_2d_array(u_old_data.group(1))
    v_old = extract_2d_array(v_old_data.group(1))
    U_old = extract_2d_array(U_old_data.group(1))
    V_old = extract_2d_array(V_old_data.group(1))
    return u, v, U, V , u_old, v_old, U_old, V_old,  p


def Recover_result(filename):
    u, v, U, V , u_old, v_old, U_old, V_old,  p = load_2d_arrays_with_numpy(filename)
    pattern = r'output_N(\d+)_Re=(\d+)_t=(\d+\.\d+)_'
    match = re.search(pattern, filename)
    if match:
        N = float(match.group(1))
        Re = float(match.group(2))
        t_0 = float(match.group(3))
        # print("N =", N)
        # print("Re =", Re)
        # print("t =", t)
    else:
        print("No match found")
    Ly = 1/N* (np.size(u,0)-2)
    Lx = 1/N* (np.size(u,1)-2)
    dy = dx = 1/N
    dt = 0.01
    data_value = np.array([Re, Ly, Lx, dy, dx ,t_0 ,dt])
    keys = ["Re", "Ly", "Lx", "dy", "dx" ,"t_0" ,"dt"]
    data_re = dict(zip(keys, data_value))
    return data_re, u, v, U, V , u_old, v_old, U_old, V_old,  p




def save_result(filename, para_name, parameter):
    with open(filename, 'w') as f:
        f.write(para_name+':\n')  # 添加注释
        f.write(','.join(map(str, parameter)) + '\n\n')  # 写入u数组,并在数组后添加换行以分隔






# def load_2d_arrays_with_numpy(filename):
#     # # filename = f'cyl_output_N=32.0, Re=150.0, t=10.0.txt'
#     # filename = f'cyl_output_N=32, Re=150, t=50.txt'
#     try:
#         with open(filename, 'r') as file:
#             content = file.read()
#     except FileNotFoundError:
#         print("文件未找到,请检查文件路径是否正确")
#         return None, None, None

#     # 辅助函数来提取并构建二维数组
#     def extract_2d_array(array_content):
#         array_data = []
#         # 匹配方括号内的所有内容,包括可能的多行
#         matches = re.findall(r'\[(.*?)]', array_content.replace('\n', ''), re.DOTALL)
#         for match in matches:
#             # 替换掉可能存在的换行符,然后解析数值
#             formatted_match = ' '.join(match.split())
#             numbers = np.fromstring(formatted_match, sep=' ')
#             array_data.append(numbers)
#         return np.array(array_data)

#     u_data = re.search(r'Array u:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
#     if not (u_data and v_data and p_data):
#         print("未能正确提取u, v, 或 p的数据。请检查数组标记和文件格式。")
#         return np.array([]), np.array([]), np.array([])

#     u = extract_2d_array(u_data.group(1))
#     v = extract_2d_array(v_data.group(1))
#     p = extract_2d_array(p_data.group(1))
#     U = extract_2d_array(U_data.group(1))
#     V = extract_2d_array(V_data.group(1))
#     u_old = extract_2d_array(u_old_data.group(1))
#     v_old = extract_2d_array(v_old_data.group(1))
#     U_old = extract_2d_array(U_old_data.group(1))
#     V_old = extract_2d_array(V_old_data.group(1))
#     return u, v, U, V , u_old, v_old, U_old, V_old,  p






def main():
    print(1)


    
    # N = 32
    # Ly = 3
    # Lx = 15

    N = 32
    Ly = 4
    Lx = 8

    dy = dx = 1/N
    Re = 150

    t = 0.01
    dt = 0.01

    global data
    data_value = np.array([Re, Ly, Lx, dy, dx ,t ,dt])
    keys = ["Re", "Ly", "Lx", "dy", "dx" ,"t" ,"dt"]
    data = dict(zip(keys, data_value))

    mesh1 =  mesher(data)
    parameter_field = mesh1.parameter_grid()

    global iF, iS
    iF = parameter_field[:,:,2]
    iS = parameter_field[:,:,3]

    outer = Circle_outer(iF, iS)
    outer_e, outer_w, outer_n, outer_s = Circle_outer_classify(outer)

    # # test ifluid
    # plt.imshow(outer_n, cmap = 'gray')
    # plt.show()



    data_re, u_0, v_0, U_0, V_0 , \
        u_old_0, v_old_0, U_old_0, V_old_0,  p \
            = Recover_result(filename = 'result/save/output_N32_Re=300_t=100.01_8x4.txt')
    
    
    drag_coeff = np.zeros(1001, dtype='float')
    lift_coeff = np.zeros_like(drag_coeff)
    for n in range(1,1001):
        time = 100+dt*n
        print(time)
        print('%.2f'%time)
        filename = 'result/save/output_N32_Re=300_t=%.2f_8x4.txt'%time
        print(f'{filename}')
        data_re, u_0, v_0, U_0, V_0 , \
        u_old_0, v_old_0, U_old_0, V_old_0,  p \
            = Recover_result(filename = 'result/save/output_N32_Re=300_t=%.2f_8x4.txt'%time)
        
        p_e = p[:,:]*outer_e[:,:]
        p_w = p[:,:]*outer_w[:,:]
        p_n = p[:,:]*outer_n[:,:]
        p_s = p[:,:]*outer_s[:,:]
        drag_coeff[n] = (np.sum(p_w) - np.sum(p_e))*dy/(1.2*0.5)
        lift_coeff[n] = (np.sum(p_s) - np.sum(p_n))*dy/(1.2*0.5)
        print(n)
    print(drag_coeff)

    


    time = np.linspace(100,110.01, 1001)

    plt.plot(time[1:], drag_coeff[1:], linestyle='-', label = 'drag coefficient')  
    # plt.plot(time, lift_coeff, linestyle='-',label = 'lift coefficient')  
    plt.title('drag coefficient')
    plt.xlabel('t')
    plt.grid(True)
    plt.legend()
    plt.savefig('drag coefficient')
    # plt.show()
    plt.close()

    plt.plot(time[1:], lift_coeff[1:], linestyle='-', label = 'lift coefficient')  
    plt.title('lift coefficient')
    plt.xlabel('t')
    plt.grid(True)
    plt.legend()
    plt.savefig('lift coefficient')
    # plt.show()
    plt.close()

    save_result('drag_coeff_8x2.txt', 'drag_coeff', drag_coeff)
    save_result('ligt_coeff_8x2.txt', 'lift_coeff', lift_coeff)

    
    

    






    



if __name__ == '__main__':
    main()

    





