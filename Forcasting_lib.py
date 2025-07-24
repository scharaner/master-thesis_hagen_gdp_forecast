import os
import pandas as pd
import numpy as np
import warnings
import urllib.request
from sklearn.preprocessing import StandardScaler


import os
import pandas as pd
import urllib.request

def load_fred_data(
    local_path: str = "./Data/FRED_QD_Data.csv",
    url: str = "https://www.stlouisfed.org/-/media/project/frbstl/stlouisfed/research/fred-qd/quarterly/current.csv"
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Lädt FRED-Quartalsdaten entweder lokal oder aus dem Internet herunter, 
    verarbeitet die Datei und gibt die Daten + Transformationscodes zurück.

    Parameter:
        local_path (str): Pfad zur lokalen CSV-Datei
        url (str): URL zur CSV-Datei im Internet

    Rückgabe:
        Tuple aus:
            - pd.DataFrame: Transformierte Daten ab Zeile 3, 'sasdate' als Index
            - pd.Series   : Transformationscodes aus der 2. Zeile
    """
    if not os.path.exists(local_path):
        print(f"Datei nicht gefunden. Lade herunter von:\n{url}\n")
        try:
            urllib.request.urlretrieve(url, local_path)
            print(f"Datei gespeichert als '{local_path}'")
        except Exception as e:
            raise RuntimeError(f"Fehler beim Herunterladen: {e}")
    else:
        print(f"Datei '{local_path}' gefunden. Lade lokal...")

    try:
        df = pd.read_csv(local_path)


        # Transformationscodes extrahieren (Zeile 2, Index = 1)
        tcodes = df.iloc[1]

        # Entferne die ersten beiden Zeilen (Metainformation + Transformationscodes)
        df = df.iloc[2:].copy()

        # Transformation: sasdate ab Zeile 3 in datetime
        df['sasdate'] = pd.to_datetime(df['sasdate'])

        # Setze Datum als Index
        df = df.set_index('sasdate')

        print("Daten erfolgreich vorbereitet.")
        #print(df.head())  # Zeige die ersten Zeilen des DataFrames an
        return df, tcodes

    except Exception as e:
        raise RuntimeError(f"Fehler beim Einlesen oder Verarbeiten der CSV-Datei: {e}")


def load_fred_MD_data(
    local_path: str = "./Data/FRED_MD_Data.csv",
    url: str = "https://www.stlouisfed.org/-/media/project/frbstl/stlouisfed/research/fred-md/quarterly/current.csv"
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Lädt FRED-Quartalsdaten entweder lokal oder aus dem Internet herunter, 
    verarbeitet die Datei und gibt die Daten + Transformationscodes zurück.

    Parameter:
        local_path (str): Pfad zur lokalen CSV-Datei
        url (str): URL zur CSV-Datei im Internet

    Rückgabe:
        Tuple aus:
            - pd.DataFrame: Transformierte Daten ab Zeile 3, 'sasdate' als Index
            - pd.Series   : Transformationscodes aus der 2. Zeile
    """
    if not os.path.exists(local_path):
        print(f"Datei nicht gefunden. Lade herunter von:\n{url}\n")
        try:
            urllib.request.urlretrieve(url, local_path)
            print(f"Datei gespeichert als '{local_path}'")
        except Exception as e:
            raise RuntimeError(f"Fehler beim Herunterladen: {e}")
    else:
        print(f"Datei '{local_path}' gefunden. Lade lokal...")

    try:
        df = pd.read_csv(local_path)


        # Transformationscodes extrahieren (Zeile 2, Index = 1)
        tcodes = df.iloc[1]

        # Entferne die ersten beiden Zeilen (Metainformation + Transformationscodes)
        df = df.iloc[2:].copy()

        # Transformation: sasdate ab Zeile 3 in datetime
        df['sasdate'] = pd.to_datetime(df['sasdate'])

        # Setze Datum als Index
        df = df.set_index('sasdate')

        print("Daten erfolgreich vorbereitet.")
        #print(df.head())  # Zeige die ersten Zeilen des DataFrames an
        return df, tcodes

    except Exception as e:
        raise RuntimeError(f"Fehler beim Einlesen oder Verarbeiten der CSV-Datei: {e}")




import numpy as np
import pandas as pd

def make_stationary(df: pd.DataFrame, tcodes: pd.Series) -> pd.DataFrame:
    """
    Transformiert ein gesamtes DataFrame gemäß eines Dictionarys von Transformationscodes 
    um die Zeitreihen stationär zu machen.

    Parameter:
        df     : pd.DataFrame – Ursprüngliches DataFrame mit Zeitreihen
        tcodes : pd.Series     – Series {Spaltenname: TCode (int)}

    Rückgabe:
        transformed_df : pd.DataFrame – Transformiertes DataFrame
    """
    transformed_cols = []
    applied_codes = {}

    for col in df.columns:
        tcode = tcodes.get(col, 1)
        s = df[col].copy()

        try:
            if tcode == 1:
                transformed = s
                desc = "Level (keine Transformation)"
            elif tcode == 2:
                transformed = s.diff()
                desc = "Erste Differenz"
            elif tcode == 3:
                transformed = s.diff().diff()
                desc = "Zweite Differenz"
            elif tcode == 4:
                transformed = np.log(s)
                desc = "Log-Level"
            elif tcode == 5:
                transformed = np.log(s).diff()
                desc = "Erste Log-Differenz (≈ Wachstumsrate)"
            elif tcode == 6:
                transformed = np.log(s).diff().diff()
                desc = "Zweite Log-Differenz"
            elif tcode == 7:
                transformed = 100 * (s / s.shift(1) - 1)
                desc = "Prozentuale Veränderung"
            else:
                transformed = s
                desc = "Unbekannter Code – keine Transformation"
        except Exception as e:
            transformed = pd.Series(np.nan, index=s.index)
            desc = f"Fehler bei Transformation: {str(e)}"

        transformed.name = col
        transformed_cols.append(transformed)
        applied_codes[col] = desc

    # Schneller und speichereffizienter Aufbau des DataFrames
    transformed_df = pd.concat(transformed_cols, axis=1)

    return transformed_df



def standardize_data(df: pd.DataFrame, columns) -> pd.DataFrame:
    """
    Standardisiert die angegebenen Spalten eines DataFrames.

    Parameter:
        df      : pd.DataFrame – DataFrame mit den zu standardisierenden Daten
        columns : list         – Liste der Spaltennamen, die standardisiert werden sollen

    Rückgabe:
        pd.DataFrame – DataFrame mit standardisierten Spalten
    """
    scaler = StandardScaler()
    df[columns] = scaler.fit_transform(df[columns])
    
    print("Standardisierung abgeschlossen.")
    #print(df[columns].head())  # Zeige die ersten Zeilen der standardisierten Spalten an
    #print(df.head())

    return df, scaler