import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pyproj import Proj

p = Proj(proj='merc', ellps='WGS84')
# начальный вектор состояния
x = np.array([0, 0, 0, 0])
# ковариционная матрица
P = np.eye(4) * 1000
# матрица эволюции
F = np.array([[1, 0, 1 / 3600, 0],
              [0, 1, 0, 1 / 3600],
              [0, 0, 1 / 3600, 0],
              [0, 0, 0, 1 / 3600]])
# матрица модели наблюдения
H = np.eye(4)
# ковариционная матрица шума процесса
Q = np.array([[2.77777778e-06, 8.33333333e-06, 1.66666667e-05, 1.66666667e-05],
              [8.33333333e-06, 2.50000000e-05, 5.00000000e-05, 5.00000000e-05],
              [1.66666667e-05, 5.00000000e-05, 1.00000000e-04, 1.00000000e-04],
              [1.66666667e-05, 5.00000000e-05, 1.00000000e-04, 1.00000000e-04]])
# ковариционная матрица шума измерения
R = np.eye(4) * 0.0001
graphs = []


def convert_coordinates(row):
    x_c, y = p(row[' Longitude'], row[' Latitude'])
    return pd.Series({'X': x_c, 'Y': y})


def load_data(path):
    df = pd.read_csv(path)
    df.loc[:232, " Bearing"] = 215.5
    df[["X", "Y"]] = df.apply(convert_coordinates, axis=1)
    df[" Bearing"] = df[" Bearing"].replace(0, np.nan)
    df[" Bearing"].interpolate(method='linear', inplace=True)
    df[" Bearing"] = np.deg2rad(df[" Bearing"])
    df['Speed (GPS)(km/h)'] = pd.to_numeric(df['Speed (GPS)(km/h)'], errors='coerce')
    df["Speed X"] = df["Speed (GPS)(km/h)"] * np.cos(df[" Bearing"])
    df["Speed Y"] = df["Speed (GPS)(km/h)"] * np.sin(df[" Bearing"])
    filtered_data = df.drop_duplicates(subset=["GPS Time"])
    scatter = go.Scatter(x=filtered_data["X"], y=filtered_data["Y"], mode='lines', line=dict(width=2, color='blue'),
                         name="Показания GPS")
    graphs.append(scatter)
    return [df, filtered_data]
    #filtered_data = df.drop_duplicates(subset=['GPS Time'])
    #iltered_data["Speed (GPS)(km/h)"] = filtered_data["Speed (GPS)(km/h)"].replace("-", np.nan)
    #iltered_data["Speed (GPS)(km/h)"].interpolate(method="linear", inplace=True)
    #iltered_data["Speed (GPS)(km/h)"] = filtered_data["Speed (GPS)(km/h)"].replace(np.nan, 13.0)
    #iltered_data["Speed (GPS)(km/h)"] = filtered_data["Speed (GPS)(km/h)"].astype(float)
    # cx = []
    # cy = []
    # R_E = 6378137
    # for i in range(len(df[" Longitude"])):
    #    cy.append(np.log(np.tan(np.pi / 4 + np.radians(df[" Latitude"][i]) / 2)) * R_E)
    #    cx.append(np.radians(df[" Latitude"][i]) * R_E)
    # df["X"] = cx
    # df["Y"] = cy

    # filtered_data = df
    #fields_gps = ["GPS Time", "X", "Y", "Speed (GPS)(km/h)", " Bearing"]
    #fields = ["GPS Time", " Latitude", " Longitude", "Speed (OBD)(km/h)", " Bearing"]
    #fields_xy = ["GPS Time", "X", "Y", "Speed (OBD)(km/h)", " Bearing"]
    #filtered_data_gps = filtered_data.loc[:, fields_gps]
    #filtered_data = filtered_data.loc[:, fields_avg]

    #filtered_data[" Bearing"] = filtered_data[" Bearing"].replace(0, np.nan)
    #filtered_data[' Bearing'].interpolate(method='linear', inplace=True)
    #filtered_data[" Bearing"] = np.deg2rad(filtered_data[" Bearing"])
    #scatter = go.Scatter(x=filtered_data["X"], y=filtered_data["Y"], mode='lines', line=dict(width=2, color='blue'),
    #                     name="Показания GPS")
    #graphs.append(scatter)
    #filtered_data_gps[" Bearing"] = filtered_data[" Bearing"]
    #filtered_data_gps["Speed X"] = filtered_data_gps["Speed (GPS)(km/h)"] * np.cos(filtered_data_gps[" Bearing"])
    #filtered_data_gps["Speed Y"] = filtered_data_gps["Speed (GPS)(km/h)"] * np.sin(filtered_data_gps[" Bearing"])
    #filtered_data["Speed X"] = filtered_data["Speed_Combined"] * np.cos(filtered_data[" Bearing"])
    #filtered_data["Speed Y"] = filtered_data["Speed_Combined"] * np.sin(filtered_data[" Bearing"])
    #return [filtered_data, filtered_data_gps]


def motion(z, x, P, is_prev):
    x = F @ x
    P = F @ P @ F.T + Q
    if not is_prev:
        x, P = update(z, x, P)
    return x, P


def update(z, x, P):
    # Шаг обновления
    y = z - H @ x
    S = H @ P @ H.T + R
    K = P @ H.T @ np.linalg.inv(S)
    x = x + K @ y
    P = P - K @ H @ P

    return x, P


vectors = []


def reset_global():
    global x
    global P
    global F
    global H
    global Q
    global R
    x = np.array([0, 0, 0, 0])
    # ковариционная матрица
    P = np.eye(4) * 1000
    # матрица эволюции
    F = np.array([[1, 0, 1 / 3600, 0],
                  [0, 1, 0, 1 / 3600],
                  [0, 0, 1 / 3600, 0],
                  [0, 0, 0, 1 / 3600]])
    # матрица модели наблюдения
    H = np.eye(4)
    # ковариционная матрица шума процесса
    Q = np.array([[2.77777778e-06, 8.33333333e-06, 1.66666667e-05, 1.66666667e-05],
                  [8.33333333e-06, 2.50000000e-05, 5.00000000e-05, 5.00000000e-05],
                  [1.66666667e-05, 5.00000000e-05, 1.00000000e-04, 1.00000000e-04],
                  [1.66666667e-05, 5.00000000e-05, 1.00000000e-04, 1.00000000e-04]])
    # ковариционная матрица шума измерения
    R = np.eye(4) * 0.0001


def sense(dataset):
    prev = dataset[0][1::]
    for z in dataset:
        global x
        global P
        if z[0] == prev[0] or "-" in z or np.nan in z:
            x, P = motion(prev, x, P, True)
        else:
            x, P = motion(z[1::], x, P, False)
        prev = z[1::]
        vectors.append(list(x))


def plot_data(**kwargs):
    layout = go.Layout(title='',
                       xaxis=dict(title=kwargs["title_x"]),
                       yaxis=dict(title=kwargs["title_y"]),
                       hovermode='closest')
    fig = go.Figure(data=kwargs["graphs"], layout=layout)
    fig.show()


if __name__ == "__main__":
    file_path = "data1.csv"
    data = load_data(file_path)
    data[0] = data[0].loc[:, ["GPS Time", "X", 'Y', 'Speed X', 'Speed Y']]
    data[1] = data[1].loc[:, ["GPS Time", "X", 'Y', 'Speed X', 'Speed Y']]
    sense(np.array(data[0]))
    X_coords = [i[0] for i in vectors]
    Y_coords = [i[1] for i in vectors]
    scatter = go.Scatter(x=X_coords, y=Y_coords, mode='lines', line=dict(width=2, color='green'), name="Оценка "
                                                                                                       "фильтра (GPS "
                                                                                                       "+ OBD)")
    graphs.append(scatter)
    vectors.clear()
    reset_global()
    sense(np.array(data[1]))
    X_coords = [i[0] for i in vectors]
    Y_coords = [i[1] for i in vectors]
    scatter = go.Scatter(x=X_coords, y=Y_coords, mode='lines', line=dict(width=2, color='red'), name="Оценка фильтра "
                                                                                                     "(GPS)")
    graphs.append(scatter)
    plot_data(title_x="X", title_y="Y", graphs=[graphs[1], graphs[2]])
    plot_data(title_x="X", title_y="Y", graphs=[graphs[0], graphs[2]])
    #plot_data(title_x="X", title_y="Y", graphs=[graphs[1], graphs[2]])
