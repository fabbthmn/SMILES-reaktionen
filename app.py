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
        # RXN unterstützt das Standard-Format Edukt>>Produkt oder Edukt>Agent>Produkt
        rxn = AllChem.ReactionFromSmarts(smiles_str, useSmiles=True)
        return AllChem.ReactionToRxnBlock(rxn)
    except:
        return None

def parse_smiles(smiles_str):
    """Unterstützt jetzt A.B >> C ODER A.B > Reagenz > C"""
    if not smiles_str or str(smiles_str).lower() == "nan": return None
    
    # Wir splitten nach ">", um zu sehen wie viele Teile wir haben
    parts = str(smiles_str).split(">")
    # Entferne leere Strings (entstehen bei >>)
    parts = [p for p in parts if p.strip() != ""]
    
    if len(parts) == 3:
        # Format: Edukt > Reagenz > Produkt
        return {"r": parts[0].split("."), "a": parts[1], "p": parts[2].split(".")}
    elif len(parts) == 2:
        # Format: Edukt >> Produkt
        return {"r": parts[0].split("."), "a": "", "p": parts[1].split(".")}
    return None

def draw_single_line(data, ratio_str, f_main, f_title, name=None):
    r_parts = str(ratio_str).replace(">>", ".").split(".")
    anzahl = len(data['r']) + len(data['p'])
    canvas_w = (anzahl * 600) + 1200 
    canvas_h = 850 # Erhöhte Höhe für mehr Platz zwischen Titel und Bild
    canvas = Image.new('RGB', (canvas_w, canvas_h), color=(255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    y_mid = canvas_h // 2
    x_cursor = 100

    # Titel (Reaktionsname)
    if name:
        draw.text((100, 30), str(name), fill="black", font=f_title)

    def draw_text_centered(txt, x, y, font):
        bbox = draw.textbbox((0, 0), txt, font=font)
        th = bbox[3] - bbox[1]
        draw.text((x, y - (th // 2) - 10), txt, fill="black", font=font)
        return bbox[2] - bbox[0]

    # --- 1. REAKTANDEN ---
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

    # --- 2. PFEIL MIT REAGENZIEN ---
    x_cursor += 60
    arrow_len = 350 # Etwas länger für Text auf dem Pfeil
    # Pfeil-Linie
    draw.line([(x_cursor, y_mid), (x_cursor + arrow_len, y_mid)], fill="black", width=4)
    # Pfeil-Spitze
    draw.polygon([(x_cursor + arrow_len, y_mid), (x_cursor + arrow_len - 30, y_mid - 15), (x_cursor + arrow_len - 30, y_mid + 15)], fill="black")
    
    # Reagenzien-Text ÜBER den Pfeil
    if data.get("a") and data["a"].strip() != "":
        f_reag = get_font(45) # Etwas kleinere Schrift als Haupt-Symbole
        reag_txt = data["a"]
        bbox_r = draw.textbbox((0, 0), reag_txt, font=f_reag)
        r_w = bbox_r[2] - bbox_r[0]
        # Text mittig über dem Pfeil platzieren
        draw.text((x_cursor + (arrow_len // 2) - (r_w // 2), y_mid - 70), reag_txt, fill="black", font=f_reag)

    x_cursor += arrow_len + 70

    # --- 3. PRODUKTE ---
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

# --- INTERFACE ---

st.title("⚗️ Chemie-Designer Pro")
t1, t2 = st.tabs(["✨ Einzel-Eingabe", "📂 Batch-Export"])

with t1:
    name = st.text_input("Reaktionsname", "Reaktion_1")
    st.info("Format-Tipp: Edukt1.Edukt2 > Reagenz > Produkt")
    c1, c2 = st.columns(2)
    s1 = c1.text_input("SMILES 1", "C*(=O)O.OCC > H2SO4 > C*(=O)OCC")
    r1 = c1.text_input("Verhältnis 1", "1>>1")
    s2 = c2.text_input("SMILES 2 (Optional)", "")
    r2 = c2.text_input("Verhältnis 2", "1>>1")

    if st.button("Generieren"):
        s_list = [s for s in [s1, s2] if s.strip()]
        img = combine_and_render(s_list, [r1, r2], name)
        if img:
            st.image(img)
            b_png = io.BytesIO(); img.save(b_png, format="PNG")
            cd1, cd2 = st.columns(2)
            cd1.download_button("💾 Bild (PNG)", b_png.getvalue(), f"{name}.png")
            rxn = export_to_rxn(s1)
            if rxn: cd2.download_button("🧬 Editierbar (RXN)", rxn, f"{name}.rxn")

with t2:
    up = st.file_uploader("Excel hochladen", type=["xlsx"])
    if up:
        df = pd.read_excel(up)
        cols = df.columns.tolist()
        cn = st.selectbox("Name Spalte", cols); cs1 = st.selectbox("SMILES A Spalte", cols)
        cs2 = st.selectbox("SMILES B Spalte (Opt)", ["Keine"] + cols); cr = st.selectbox("Verhältnis Spalte", ["Keine"] + cols)

        if st.button("🚀 Reaktionen generieren"):
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                for idx, row in df.iterrows():
                    nom_raw = str(row[cn])
                    nom = nom_raw.replace("/", "_").replace("\\", "_")
                    s_list = [str(row[cs1])]
                    ratio = str(row[cr]) if cr != "Keine" else "1>>1"
                    if cs2 != "Keine" and pd.notna(row[cs2]): s_list.append(str(row[cs2]))
                    
                    img = combine_and_render(s_list, [ratio, ratio], nom_raw)
                    if img:
                        st.markdown(f"### {idx+1}. {nom_raw}")
                        st.image(img)
                        img_byte_arr = io.BytesIO(); img.save(img_byte_arr, format='PNG')
                        rxn_data = export_to_rxn(str(row[cs1]))
                        c_d1, c_d2 = st.columns(2)
                        c_d1.download_button(f"PNG: {nom}", img_byte_arr.getvalue(), f"{nom}.png", key=f"p_{idx}")
                        if rxn_data:
                            c_d2.download_button(f"RXN: {nom}", rxn_data, f"{nom}.rxn", key=f"r_{idx}")
                        zip_file.writestr(f"{nom}.png", img_byte_arr.getvalue())
                        if rxn_data: zip_file.writestr(f"{nom}.rxn", rxn_data)
                        st.divider()
            
            st.success("✅ Alle Reaktionen fertig erstellt!")
            st.download_button(
                label="📦 ALLE REAKTIONEN ALS ZIP HERUNTERLADEN",
                data=zip_buffer.getvalue(),
                file_name="chemie_batch_export.zip",
                mime="application/zip",
                use_container_width=True
            )