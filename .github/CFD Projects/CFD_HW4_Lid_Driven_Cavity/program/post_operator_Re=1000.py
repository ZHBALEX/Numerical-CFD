import numpy as np
import re
import matplotlib.pyplot as plt

def load_2d_arrays_with_numpy(label):
    import os
    import sys

    # 将工作目录改变到脚本所在目录
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    filename = f'data/Re1000_output_N={label}.txt'
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

    # 改进正则表达式以匹配整个数组块
    u_data = re.search(r'Array u:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    v_data = re.search(r'Array v:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    p_data = re.search(r'Array p:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)

    if not (u_data and v_data and p_data):
        print("未能正确提取u, v, 或 p的数据。请检查数组标记和文件格式。")
        return np.array([]), np.array([]), np.array([])

    u = extract_2d_array(u_data.group(1))
    v = extract_2d_array(v_data.group(1))
    p = extract_2d_array(p_data.group(1))

    return u, v, p

def smalluvp(u,v,p):
    return u[1:-1,1:-1],v[1:-1,1:-1],p[1:-1,1:-1] 

def Center(u,v,p, label):
    u,v,p = smalluvp(u,v,p)
    
    center = int(label/2)
    u_x_center = (u[:,center] + u[:,center-1])/2
    v_x_center = (v[:,center] + v[:,center-1])/2
    p_x_center = (p[:,center] + p[:,center-1])/2

    u_y_center = (u[center,:] + u[center-1,:])/2
    v_y_center = (v[center,:] + v[center-1,:])/2
    p_y_center = (p[center,:] + p[center-1,:])/2

    return u_x_center, v_x_center, p_x_center\
    ,u_y_center, v_y_center, p_y_center


# def drawing(u_x_c, u_y_c):
    
    
def get_ue_x_c():
    ye_0 = [1.00000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]
    ue_x_c_0 = [1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.06080, -0.10648, -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0.00000]
    return ye_0[::-1], ue_x_c_0[::-1]

def get_ue_y_c():
    xe_0 = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000]
    xe_y_c_0 = [0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, -0.42665, -0.31966, 0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.00000]
    return xe_0[::-1], xe_y_c_0[::-1]


# def drawing(u_x_c_16, u_x_c_64, u_x_c_128):
#     ye, ue_x_c = get_ue_x_c()
#     y = np.linspace(0,N*1/N, N)
#     plt.figure(figsize=(8, 8))
#     plt.plot(ye, ue_x_c, label='Paper Result')
#     plt.plot(y, u_x_c, label='Result')
#     plt.xlabel('X')
#     plt.ylabel('u velocity')
#     plt.legend()
#     # plt.title('u Contour')
#     # plt.savefig(f'u_contour_{label}.png')
#     plt.show()


def drawing_x_c(u_x_c_16, u_x_c_64, u_x_c_128):
    ye, ue_x_c_vertical = get_ue_x_c()
    y_16 = np.linspace(0,1, 16)
    y_64 = np.linspace(0,1, 64)
    y_128 = np.linspace(0,1, 128)
    plt.figure(figsize=(8, 8))
    plt.plot( ue_x_c_vertical, ye, 'X', label='Reference', color= '#0072B2')
    plt.plot( u_x_c_16,y_16, '*', label='N_16',color= '#56B4E9')
    plt.plot( u_x_c_64,y_64, '+', label='N_64',color= '#D55E00')
    plt.plot( u_x_c_128,y_128, '-', label='N_128',color= '#CC79A7')
    plt.xlabel('Y')
    plt.ylabel('u velocity')
    plt.legend()
    plt.title('Re=1000_u along verticle line through geometric center')
    plt.savefig('Re=1000_u along verticle line through geometric center')
    # plt.show()



    

def drawing_y_c(v_y_c_16, v_y_c_64, v_y_c_128):
    xe, ue_y_c_vertical = get_ue_y_c()
    plt.figure(figsize=(12,12))
    x_16 = np.linspace(0,16*1/16, 16)
    x_64 = np.linspace(0,64*1/64, 64)
    x_128 = np.linspace(0,128*1/128, 128)

    plt.plot(xe, ue_y_c_vertical, 'X', label='Reference',color= '#0072B2')
    plt.plot(x_16, v_y_c_16, '*', label='N_16',color= '#56B4E9')
    plt.plot(x_64, v_y_c_64, '+', label='N_64',color= '#D55E00')
    plt.plot(x_128, v_y_c_128, '-', label='N_128',color= '#CC79A7')
    plt.xlabel('X')
    plt.ylabel('v velocity')
    plt.legend()
    plt.title('Re=1000_v along horizontal line through geometric center')
    plt.savefig('Re=1000_v along horizontal line through geometric center')
    # plt.show()


def center_get(N):
    u, v, p = load_2d_arrays_with_numpy(N)
    return Center(u,v,p, N)




def draw_streamline(u, v, len_x, len_y ,label):
    rows, cols = np.shape(u)  
    X = np.linspace(0, len_x, cols)
    Y = np.linspace(0, len_y, rows)
    plt.figure(figsize=(12, 12))
    plt.axis([0,1,0,1])
    plt.streamplot(X, Y, u, v, density= 3)
    plt.savefig(f'Re1000_streamline_{label}.png')
    # plt.show()



def main():
    # 使用示例
    # N = 128
    # u, v, p = load_2d_arrays_with_numpy(N)
    # print(u[-2,:])
    # u_x_c, v_x_c, p_x_c, u_y_c, v_y_c, p_y_c = Center(u,v,p, N)

    u_x_c_128, v_x_c_128, p_x_c_128, u_y_c_128, v_y_c_128, p_y_c_128 = center_get(128)
    u_x_c_64, v_x_c_64, p_x_c_64, u_y_c_64, v_y_c_64, p_y_c_64 = center_get(64)
    u_x_c_16, v_x_c_16, p_x_c_16, u_y_c_16, v_y_c_16, p_y_c_16 = center_get(16)
    
    # N=128
    # u, v, p = load_2d_arrays_with_numpy(N)
    # draw_streamline(u[1:-1,1:-1],v[1:-1,1:-1], 1, 1 , N)

    drawing_x_c(u_x_c_16, u_x_c_64, u_x_c_128)
    drawing_y_c(v_y_c_16, v_y_c_64, v_y_c_128)



if __name__ == '__main__':
    main()

    






