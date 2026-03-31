import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io
import zipfile

# --- 1. GLOBALE KONFIGURATION ---
st.set_page_config(page_title="Chemie-Designer Pro", page_icon="🧪", layout="wide")

# --- 2. HILFSFUNKTIONEN ---

def get_font(size):
    font_paths = [
        "arial.ttf", 
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 
        "/System/Library/Fonts/Helvetica.ttc"
    ]
    for path in font_paths:
        try:
            return ImageFont.truetype(path, size)
        except:
            continue
    return ImageFont.load_default()

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
    opts.bondLineWidth = 3
    opts.minFontSize = 25
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))

def export_to_rxn(smiles_str):
    try:
        rxn = AllChem.ReactionFromSmarts(smiles_str, useSmiles=True)
        return AllChem.ReactionToRxnBlock(rxn)
    except:
        return None

def parse_smiles(smiles_str):
    if not smiles_str or str(smiles_str).lower() == "nan": return None
    if ">>" in smiles_str:
        parts = smiles_str.split(">>")
        return {"r": parts[0].split("."), "p": parts[1].split(".")}
    return None

def draw_single_line(data, ratio_str, f_main, f_title, name=None):
    r_parts = str(ratio_str).replace(">>", ".").split(".")
    anzahl = len(data['r']) + len(data['p'])
    canvas_w = (anzahl * 600) + 1200 
    canvas_h = 700
    canvas = Image.new('RGB', (canvas_w, canvas_h), color=(255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    y_mid = canvas_h // 2
    x_cursor = 100

    if name:
        draw.text((100, 40), str(name), fill="black", font=f_title)

    def draw_text_centered(txt, x, y, font):
        bbox = draw.textbbox((0, 0), txt, font=font)
        th = bbox[3] - bbox[1]
        draw.text((x, y - (th // 2) - 10), txt, fill="black", font=font)
        return bbox[2] - bbox[0]

    for i, sm in enumerate(data['r']):
        coeff = r_parts[i].strip() if i < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            x_cursor += draw_text_centered(coeff, x_cursor, y_mid, f_main) + 30
        m_img = generate_mol_img(sm)
        if m_img:
            canvas.paste(m_img, (x_cursor, y_mid - 250), m_img.convert("RGBA"))
            x_cursor += 500
        if i < len(data['r']) - 1:
            x_cursor += 30
            x_cursor += draw_text_centered("+", x_cursor, y_mid, f_main) + 50

    x_cursor += 60
    draw.line([(x_cursor, y_mid), (x_cursor + 250, y_mid)], fill="black", width=4)
    draw.polygon([(x_cursor + 250, y_mid), (x_cursor + 220, y_mid - 15), (x_cursor + 220, y_mid + 15)], fill="black")
    x_cursor += 320

    offset = len(data['r'])
    for i, sm in enumerate(data['p']):
        c_idx = offset + i
        coeff = r_parts[c_idx].strip() if c_idx < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            x_cursor += draw_text_centered(coeff, x_cursor, y_mid, f_main) + 30
        m_img = generate_mol_img(sm)
        if m_img:
            canvas.paste(m_img, (x_cursor, y_mid - 250), m_img.convert("RGBA"))
            x_cursor += 500
        if i < len(data['p']) - 1:
            x_cursor += 30
            x_cursor += draw_text_centered("+", x_cursor, y_mid, f_main) + 50

    return canvas.crop((0, 0, x_cursor + 150, canvas_h))

def combine_and_render(smiles_list, ratio_list, reaction_name):
    f_title = get_font(60)
    f_main = get_font(90)
    images = []
    for i, sm in enumerate(smiles_list):
        data = parse_smiles(sm)
        if data:
            images.append(draw_single_line(data, ratio_list[i], f_main, f_title, reaction_name if i == 0 else None))
    if not images: return None
    max_w = max(img.width for img in images)
    total_h = sum(img.height for img in images) + 20
    combined = Image.new('RGB', (max_w, total_h), (255, 255, 255))
    y_off = 0
    for img in images:
        combined.paste(img, (0, y_off))
        y_off += img.height
    return combined

# --- 4. INTERFACE ---

st.title("⚗️ Chemie-Designer Pro")
t1, t2 = st.tabs(["✨ Einzel-Eingabe", "📂 Batch-Export"])

with t1:
    name = st.text_input("Reaktionsname", "Reaktion_1")
    c1, c2 = st.columns(2)
    s1 = c1.text_input("SMILES 1", "C*(=O)O.OCC>>C*(=O)OCC")
    r1 = c1.text_input("Verhältnis 1", "1>>1")
    s2 = c2.text_input("SMILES 2 (Optional)", "")
    r2 = c2.text_input("Verhältnis 2", "1>>1")

    if st.button("Generieren"):
        s_list = [s for s in [s1, s2] if s.strip()]
        img = combine_and_render(s_list, [r1, r2], name)
        if img:
            st.image(img)
            b_png = io.BytesIO()
            img.save(b_png, format="PNG")
            col_d1, col_d2 = st.columns(2)
            col_d1.download_button("💾 Bild (PNG)", b_png.getvalue(), f"{name}.png")
            rxn_data = export_to_rxn(s1)
            if rxn_data:
                col_d2.download_button("🧬 Editierbar (RXN)", rxn_data, f"{name}.rxn")

with t2:
    up = st.file_uploader("Excel hochladen", type=["xlsx"])
    if up:
        df = pd.read_excel(up)
        cols = df.columns.tolist()
        cn = st.selectbox("Name", cols); cs1 = st.selectbox("SMILES A", cols)
        cs2 = st.selectbox("SMILES B (Opt)", ["Keine"] + cols); cr = st.selectbox("Verhältnis", ["Keine"] + cols)

        if st.button("🚀 Alle Generieren & ZIP vorbereiten"):
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED, False) as zip_file:
                for idx, row in df.iterrows():
                    nom = str(row[cn]).replace("/", "_") # Dateinamen sicher machen
                    s_list = [str(row[cs1])]
                    ratio = str(row[cr]) if cr != "Keine" else "1>>1"
                    if cs2 != "Keine" and pd.notna(row[cs2]): s_list.append(str(row[cs2]))
                    
                    img = combine_and_render(s_list, [ratio, ratio], nom)
                    if img:
                        st.write(f"Verarbeite: {nom}")
                        # Bild für ZIP
                        img_byte_arr = io.BytesIO()
                        img.save(img_byte_arr, format='PNG')
                        zip_file.writestr(f"{nom}.png", img_byte_arr.getvalue())
                        
                        # RXN für ZIP
                        rxn_data = export_to_rxn(str(row[cs1]))
                        if rxn_data:
                            zip_file.writestr(f"{nom}.rxn", rxn_data)
            
            st.success("✅ Alle Reaktionen verarbeitet!")
            st.download_button(
                label="📦 ALLES HERUNTERLADEN (ZIP)",
                data=zip_buffer.getvalue(),
                file_name="chemie_export.zip",
                mime="application/zip",
                use_container_width=True
            )