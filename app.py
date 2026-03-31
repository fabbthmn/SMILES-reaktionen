import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io

# --- 1. GLOBALE KONFIGURATION ---
st.set_page_config(page_title="Chemie-Designer Pro", page_icon="🧪", layout="wide")

# --- 2. FUNKTIONEN ---
def generate_mol_img(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    
    # R-Label Logik
    dummy_idx = 0
    labels = ["R", "R'", "R''", "R'''"]
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "*":
            label = labels[dummy_idx] if dummy_idx < len(labels) else f"R{dummy_idx}"
            atom.SetProp("atomLabel", label)
            dummy_idx += 1

    AllChem.Compute2DCoords(mol)
    
    # GRÖSSERES BILD: 500x500 statt 400x400
    d = rdMolDraw2D.MolDraw2DCairo(500, 500)
    opts = d.drawOptions()
    opts.bondLineWidth = 4         # Noch dickere Bindungen für die Optik
    opts.minFontSize = 25          # Atom-Labels (O, N, R) werden größer
    opts.fixedBondLength = 45      # Längere Bindungen lassen das Molekül größer wirken
    
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))

def export_to_rxn(smiles_str):
    """Erstellt einen MDL RXN Block für KingDraw/BioVIA."""
    try:
        rxn = AllChem.ReactionFromSmarts(smiles_str, useSmiles=True)
        return AllChem.ReactionToRxnBlock(rxn)
    except:
        return None

def parse_smiles(smiles_str):
    if not smiles_str or str(smiles_str).lower() == "nan": return None
    if ">>" in smiles_str:
        parts = smiles_str.split(">>")
        reactants = parts[0].split(".")
        products = parts[1].split(".")
        # Agent-Check (falls Format A.B>Agent>C)
        return {"r": reactants, "a": "", "p": products}
    return None

def draw_reaction_line(data, ratio_str, reaction_name="Unbenannte Reaktion"):
    """Zeichnet die Reaktion mit einer dynamisch berechneten Breite."""
    r_parts = str(ratio_str).replace(">>", ".").split(".")
    
    # --- DYNAMISCHE BREITE BERECHNEN ---
    # Wir schätzen: pro Molekül ca. 550px, plus Pfeil (400px), plus Puffer
    anzahl_elemente = len(data['r']) + len(data['p'])
    dynamische_breite = (anzahl_elemente * 550) + 800 
    # Sicherstellen, dass die Leinwand nicht zu klein für den Namen ist
    dynamische_breite = max(dynamische_breite, 1500) 
    
    canvas_h = 700
    canvas = Image.new('RGB', (dynamische_breite, canvas_h), color=(255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    
    # --- SCHRIFTARTEN (wie besprochen) ---
    font_title = ImageFont.load_default()
    font_main = ImageFont.load_default()
    try: 
        font_title = ImageFont.truetype("arial.ttf", 45)
        font_main = ImageFont.truetype("arial.ttf", 65)
    except: 
        try:
            font_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 45)
            font_main = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 65)
        except:
            pass

    draw.text((50, 30), str(reaction_name), fill="black", font=font_title)
    
    x_cursor = 50
    y_midline = canvas_h // 2  # Automatisch 350

    def draw_centered_text(text, x, y_mid, font, color="black"):
        bbox = draw.textbbox((0, 0), text, font=font)
        text_h = bbox[3] - bbox[1]
        draw.text((x, y_mid - (text_h // 2) - 10), text, fill=color, font=font)
        return bbox[2] - bbox[0] # Return width

    # 1. REAKTANDEN
    for i, sm in enumerate(data['r']):
        coeff = r_parts[i].strip() if i < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            w = draw_centered_text(coeff, x_cursor, y_midline, font_main)
            x_cursor += w + 20
        
        mol_img = generate_mol_img(sm)
        if mol_img:
            mol_y = y_midline - (mol_img.height // 2)
            canvas.paste(mol_img, (x_cursor, mol_y), mol_img.convert("RGBA"))
            x_cursor += mol_img.width
        
        if i < len(data['r']) - 1:
            x_cursor += 20
            w = draw_centered_text("+", x_cursor, y_midline, font_main)
            x_cursor += w + 40

    # 2. PFEIL
    x_cursor += 40
    arrow_y = y_midline
    draw.line([(x_cursor, arrow_y), (x_cursor + 250, arrow_y)], fill="black", width=5)
    draw.polygon([(x_cursor + 250, arrow_y), (x_cursor + 230, arrow_y - 15), (x_cursor + 230, arrow_y + 15)], fill="black")
    x_cursor += 300

    # 3. PRODUKTE
    offset = len(data['r'])
    for i, sm in enumerate(data['p']):
        coeff_idx = offset + i
        coeff = r_parts[coeff_idx].strip() if coeff_idx < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            w = draw_centered_text(coeff, x_cursor, y_midline, font_main)
            x_cursor += w + 20
        
        mol_img = generate_mol_img(sm)
        if mol_img:
            mol_y = y_midline - (mol_img.height // 2)
            canvas.paste(mol_img, (x_cursor, mol_y), mol_img.convert("RGBA"))
            x_cursor += mol_img.width
        
        if i < len(data['p']) - 1:
            x_cursor += 20
            w = draw_centered_text("+", x_cursor, y_midline, font_main)
            x_cursor += w + 40

    return canvas.crop((0, 0, x_cursor + 50, canvas_h))

# --- 3. INTERFACE ---

st.title("⚗️ Chemie-Designer Pro")
tab1, tab2 = st.tabs(["✨ Einzel-Eingabe", "📂 Batch-Export"])

with tab1:
    col1, col2 = st.columns([2, 1])
    with col1:
        name = st.text_input("Reaktionsname", "Veresterung")
        smiles = st.text_input("SMILES (A.B>>C)", "C*(=O)O.OCC>>C*(=O)OCC")
        ratio = st.text_input("Verhältnisse (z.B. 2>>1)", "1>>1")

    if st.button("Generieren"):
        data = parse_smiles(smiles)
        if data:
            img = draw_reaction_line(data, ratio, name)
            st.image(img)
            
            # Downloads
            c1, c2 = st.columns(2)
            # Bild
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            c1.download_button("💾 Bild speichern (PNG)", buf.getvalue(), f"{name}.png")
            
            # Chemie-Datei für KingDraw
            rxn_data = export_to_rxn(smiles)
            if rxn_data:
                c2.download_button("🧬 Editierbar (KingDraw/BioVIA)", rxn_data, f"{name}.rxn", "chemical/x-mdl-rxnfile")
        else:
            st.error("SMILES Format ungültig. Nutze '>>' als Trenner.")

# --- TAB 2: EXCEL-VERARBEITUNG ---
with tab2:
    uploaded_file = st.file_uploader("Excel-Datei hochladen (.xlsx)", type=["xlsx"])
    
    if uploaded_file:
        df = pd.read_excel(uploaded_file)
        st.write("### 1. Daten-Vorschau")
        st.dataframe(df.head(3))
        
        st.divider()
        st.write("### 2. Spalten-Zuordnung")
        col_a, col_b, col_c = st.columns(3)
        
        # Hier wählst du aus, welche Spalte was ist
        all_cols = df.columns.tolist()
        col_smiles = col_a.selectbox("Spalte mit SMILES (A.B>>C)", all_cols)
        col_name = col_b.selectbox("Spalte mit Namen", all_cols)
        col_ratio = col_c.selectbox("Spalte mit Verhältnissen (optional)", ["Keine"] + all_cols)

        if st.button("🚀 Alle Reaktionen generieren", key="btn_batch"):
            for idx, row in df.iterrows():
                # Daten aus der gewählten Spalte ziehen
                nom = str(row[col_name])
                smiles_val = str(row[col_smiles])
                
                # Verhältnis prüfen
                if col_ratio != "Keine":
                    ratio_val = str(row[col_ratio])
                else:
                    ratio_val = "1>>1"
                
                data = parse_smiles(smiles_val)
                
                if data:
                    st.write(f"---")
                    st.subheader(f"Reaktion: {nom}")
                    
                    # Bild generieren
                    img = draw_reaction_line(data, ratio_val, nom)
                    st.image(img)
                    
                    # Download Buttons
                    c1, c2 = st.columns(2)
                    
                    buf = io.BytesIO()
                    img.save(buf, format="PNG")
                    c1.download_button(f"PNG: {nom}", buf.getvalue(), f"{nom}.png", key=f"p_{idx}")
                    
                    rxn_data = export_to_rxn(smiles_val)
                    if rxn_data:
                        c2.download_button(f"KingDraw: {nom}", rxn_data, f"{nom}.rxn", key=f"k_{idx}")
                else:
                    st.warning(f"⚠️ Zeile {idx}: SMILES '{smiles_val}' konnte nicht gelesen werden.")