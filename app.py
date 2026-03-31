import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io
import os

# --- 1. GLOBALE KONFIGURATION ---
st.set_page_config(page_title="Chemie-Designer Pro", page_icon="🧪", layout="wide")

# --- 2. HILFSFUNKTION FÜR SCHRIFTEN ---
def get_huge_font(size):
    """Versucht eine skalierbare Schriftart zu laden, da Default-Fonts nicht skalierbar sind."""
    font_paths = [
        "arial.ttf",                          # Windows/Lokal
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", # Linux/Streamlit Cloud
        "/System/Library/Fonts/Helvetica.ttc", # Mac
        "LiberationSans-Regular.ttf"           # Alternative Linux
    ]
    for path in font_paths:
        try:
            return ImageFont.truetype(path, size)
        except:
            continue
    return ImageFont.load_default() # Letzter Ausweg (leider nicht skalierbar)

# --- 3. KERN-FUNKTIONEN ---

def generate_mol_img(smiles):
    if not smiles or str(smiles).lower() == "nan": return None
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    
    dummy_idx = 0
    labels = ["R", "R'", "R''", "R'''"]
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "*":
            label = labels[dummy_idx] if dummy_idx < len(labels) else f"R{dummy_idx}"
            atom.SetProp("atomLabel", label)
            dummy_idx += 1

    AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(500, 500) 
    opts = d.drawOptions()
    opts.bondLineWidth = 4
    opts.minFontSize = 30
    opts.fixedBondLength = 45
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))

def parse_smiles(smiles_str):
    if not smiles_str or str(smiles_str).lower() == "nan": return None
    if ">>" in smiles_str:
        parts = smiles_str.split(">>")
        return {"r": parts[0].split("."), "p": parts[1].split(".")}
    return None

def draw_single_reaction_line(data, ratio_str, f_main, f_title, reaction_name=None):
    r_parts = str(ratio_str).replace(">>", ".").split(".")
    
    anzahl_elemente = len(data['r']) + len(data['p'])
    canvas_w = (anzahl_elemente * 600) + 1000
    canvas_h = 800
    
    canvas = Image.new('RGB', (canvas_w, canvas_h), color=(255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    y_mid = canvas_h // 2
    x_cursor = 100

    # TITEL (Falls vorhanden)
    if reaction_name:
        draw.text((100, 50), str(reaction_name), fill="black", font=f_title)

    def draw_val(txt, x, y, font):
        bbox = draw.textbbox((0, 0), txt, font=font)
        tw = bbox[2] - bbox[0]
        th = bbox[3] - bbox[1]
        draw.text((x, y - (th // 2) - 15), txt, fill="black", font=font)
        return tw

    # REAKTANDEN
    for i, sm in enumerate(data['r']):
        coeff = r_parts[i].strip() if i < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            x_cursor += draw_val(coeff, x_cursor, y_mid, f_main) + 30
        
        m_img = generate_mol_img(sm)
        if m_img:
            canvas.paste(m_img, (x_cursor, y_mid - 250), m_img.convert("RGBA"))
            x_cursor += 500
        
        if i < len(data['r']) - 1:
            x_cursor += 20
            x_cursor += draw_val("+", x_cursor, y_mid, f_main) + 40

    # PFEIL (Wieder normaler)
    x_cursor += 50
    draw.line([(x_cursor, y_mid), (x_cursor + 250, y_mid)], fill="black", width=5)
    draw.polygon([(x_cursor + 250, y_mid), (x_cursor + 225, y_mid - 15), (x_cursor + 225, y_mid + 15)], fill="black")
    x_cursor += 300

    # PRODUKTE
    offset = len(data['r'])
    for i, sm in enumerate(data['p']):
        c_idx = offset + i
        coeff = r_parts[c_idx].strip() if c_idx < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            x_cursor += draw_val(coeff, x_cursor, y_mid, f_main) + 30
        
        m_img = generate_mol_img(sm)
        if m_img:
            canvas.paste(m_img, (x_cursor, y_mid - 250), m_img.convert("RGBA"))
            x_cursor += 500
        
        if i < len(data['p']) - 1:
            x_cursor += 20
            x_cursor += draw_val("+", x_cursor, y_mid, f_main) + 40

    return canvas.crop((0, 0, x_cursor + 100, canvas_h))

def combine_reactions(smiles_list, ratio_list, reaction_name):
    # SCHRIFTEN EXTREM GROSS DEFINIEREN
    f_title = get_huge_font(90)
    f_main = get_huge_font(120)

    images = []
    for i, sm in enumerate(smiles_list):
        data = parse_smiles(sm)
        if data:
            images.append(draw_single_reaction_line(data, ratio_list[i], f_main, f_title, reaction_name if i == 0 else None))
    
    if not images: return None

    max_w = max(img.width for img in images)
    total_h = sum(img.height for img in images)
    combined = Image.new('RGB', (max_w, total_h), (255, 255, 255))
    
    y_off = 0
    for img in images:
        combined.paste(img, (0, y_off))
        y_off += img.height
    return combined

# --- 4. INTERFACE ---
st.title("⚗️ Chemie-Designer Pro")
tab1, tab2 = st.tabs(["✨ Einzel-Eingabe", "📂 Batch-Export"])

with tab1:
    name = st.text_input("Reaktionsname", "Veresterung")
    c1, c2 = st.columns(2)
    sm1 = c1.text_input("SMILES 1", "C*(=O)O.OCC>>C*(=O)OCC")
    ra1 = c1.text_input("Verhältnis 1", "1>>1")
    sm2 = c2.text_input("SMILES 2 (Optional)", "")
    ra2 = c2.text_input("Verhältnis 2", "1>>1")

    if st.button("Generieren"):
        s_list = [s for s in [sm1, sm2] if s.strip()]
        img = combine_reactions(s_list, [ra1, ra2], name)
        if img:
            st.image(img)
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            st.download_button("💾 Download", buf.getvalue(), f"{name}.png")

with tab2:
    up = st.file_uploader("Excel hochladen", type=["xlsx"])
    if up:
        df = pd.read_excel(up)
        cols = df.columns.tolist()
        c_n = st.selectbox("Name Spalte", cols)
        c_s1 = st.selectbox("SMILES 1 Spalte", cols)
        c_s2 = st.selectbox("SMILES 2 Spalte (Optional)", ["Keine"] + cols)
        c_r = st.selectbox("Verhältnis Spalte", ["Keine"] + cols)

        if st.button("🚀 Batch Start"):
            for idx, row in df.iterrows():
                nom = str(row[c_n])
                s_list = [str(row[c_s1])]
                ratio = str(row[c_r]) if c_r != "Keine" else "1>>1"
                if c_s2 != "Keine" and pd.notna(row[c_s2]):
                    s_list.append(str(row[c_s2]))
                
                img = combine_reactions(s_list, [ratio, ratio], nom)
                if img:
                    st.image(img)
                    buf = io.BytesIO()
                    img.save(buf, format="PNG")
                    st.download_button(f"Save {nom}", buf.getvalue(), f"{nom}.png", key=f"b_{idx}")