import streamlit as st
import pandas as pd
import requests

# --- API 함수 (기존과 동일, 안정성을 위해 timeout 유지) ---
def get_pharos_data(gene_symbols):
    url = 'https://pharos-api.ncats.io/graphql'
    results = {}
    query = """
    query getTarget($gene: String!) {
      target(q: { sym: $gene }) {
        sym, name, tdl, fam, uniprot
        publications(top: 10) { pmid, title, journal, date }
      }
    }
    """
    for gene in gene_symbols:
        try:
            response = requests.post(url, json={'query': query, 'variables': {'gene': gene}}, timeout=15)
            if response.status_code == 200:
                data = response.json().get('data', {}).get('target')
                results[gene] = data if data else {'error': 'No data'}
            else:
                results[gene] = {'error': f'HTTP {response.status_code}'}
        except:
            results[gene] = {'error': 'Connection Error'}
    return results

def get_opentargets_data(uniprot_id):
    if not uniprot_id: return None
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    query = """
    query targetByUniprot($uId: [String!]!) {
      mapIds(queryTerms: $uId) {
        mappings { hits { object { ... on Target {
          knownDrugs { count, rows { drug { name }, phase, status } }
          associatedDiseases { count }
        } } } }
      }
    }
    """
    try:
        res = requests.post(url, json={'query': query, 'variables': {"uId": [uniprot_id]}}, timeout=15)
        return res.json()['data']['mapIds']['mappings'][0]['hits'][0]['object']
    except: return None

def get_alphafold_pdb(uniprot_id):
    if not uniprot_id: return None, None
    try:
        res = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}", timeout=15)
        pdb_url = res.json()[0]['pdbUrl']
        return pdb_url, requests.get(pdb_url).content
    except: return None, None

# --- UI 설정 ---
st.set_page_config(page_title="Biobytes Analyzer", layout="wide")
st.title("🧬 Biobytes Target Analyzer")

# 1. 입력 영역 (고정)
input_text = st.text_input("유전자 기호 입력 (쉼표 구분)", key="main_input")
analyze_button = st.button("데이터 분석 및 PDB 찾기")

# 2. 결과 표시를 위한 미리 정의된 컨테이너 영역 (핵심!)
result_area = st.container()

if analyze_button and input_text:
    gene_list = [g.strip().upper() for g in input_text.split(",") if g.strip()]
    
    with st.spinner('데이터 수집 중...'):
        pharos_results = get_pharos_data(gene_list)
        
        # 결과 영역 안에서 순차적으로 그리기
        with result_area:
            for gene in gene_list:
                info = pharos_results.get(gene, {})
                
                if not info or 'error' in info:
                    st.warning(f"⚠️ {gene}: 데이터를 찾을 수 없습니다.")
                    continue
                
                uniprot_id = info.get('uniprot')
                ot_data = get_opentargets_data(uniprot_id)
                pdb_url, pdb_content = get_alphafold_pdb(uniprot_id)
                
                # --- 개별 카드 렌더링 ---
                st.markdown(f"## 🎯 {gene} Report")
                
                # 메트릭 섹션
                c1, c2, c3 = st.columns(3)
                c1.metric("Family", info.get('fam', 'N/A'))
                c2.metric("TDL", info.get('tdl', 'N/A'))
                c3.metric("Drugs", ot_data['knownDrugs']['count'] if ot_data else 0)

                # 정보 그리드 (Expander 제거하여 델타 인덱스 에러 방지)
                col_left, col_right = st.columns(2)
                
                with col_left:
                    st.write("**💊 임상 약물 Top 5**")
                    if ot_data and ot_data['knownDrugs']['count'] > 0:
                        rows = ot_data['knownDrugs']['rows'][:5]
                        drug_list = [{"Name": r['drug']['name'], "Phase": r['phase']} for r in rows]
                        st.table(pd.DataFrame(drug_list))
                    else:
                        st.write("정보 없음")

                with col_right:
                    st.write("**🔬 AlphaFold PDB**")
                    if pdb_url and pdb_content:
                        st.download_button(
                            label=f"💾 {gene} PDB Download",
                            data=pdb_content,
                            file_name=f"{gene}.pdb",
                            key=f"btn_{gene}" # 단순하고 명확한 키
                        )
                    else:
                        st.write("구조 데이터 없음")
                
                st.markdown("---")
