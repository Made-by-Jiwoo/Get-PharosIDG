import streamlit as st
import pandas as pd
import requests

# --- Pharos API 호출 함수 ---
def get_pharos_data(gene_symbols):
    url = 'https://pharos-api.ncats.io/graphql'
    results = {}

    query = """
    query getTarget($gene: String!) {
      target(q: { sym: $gene }) {
        sym
        name
        tdl
        fam
        uniprot
        publications(top: 10) {
          pmid
          title
          journal
          date
        }
      }
    }
    """

    for gene in gene_symbols:
        variables = {'gene': gene}
        try:
            response = requests.post(url, json={'query': query, 'variables': variables}, timeout=15)
            if response.status_code == 200:
                resp_json = response.json()
                if 'errors' in resp_json:
                    results[gene] = {'error': resp_json['errors'][0]['message']}
                    continue
                
                target_data = resp_json.get('data', {}).get('target')
                if target_data:
                    results[gene] = target_data
                else:
                    results[gene] = {'error': 'No data found in Pharos'}
            else:
                results[gene] = {'error': f'HTTP {response.status_code}'}
        except Exception as e:
            results[gene] = {'error': str(e)}
            
    return results

# --- Open Targets API 호출 함수 ---
def get_opentargets_data(uniprot_id):
    if not uniprot_id:
        return None
    
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    query = """
    query targetByUniprot($uId: [String!]!) {
      mapIds(queryTerms: $uId) {
        mappings {
          hits {
            object {
              ... on Target {
                id
                approvedSymbol
                knownDrugs {
                  count
                  rows {
                    drug { name }
                    phase
                    status
                  }
                }
                associatedDiseases {
                  count
                }
              }
            }
          }
        }
      }
    }
    """
    
    variables = {"uId": [uniprot_id]}
    try:
        response = requests.post(url, json={'query': query, 'variables': variables}, timeout=15)
        if response.status_code == 200:
            res_json = response.json()
            mappings = res_json.get('data', {}).get('mapIds', {}).get('mappings', [])
            if mappings and mappings[0].get('hits'):
                return mappings[0]['hits'][0].get('object')
    except:
        pass
    return None

# --- AlphaFold API 호출 함수 ---
def get_alphafold_pdb(uniprot_id):
    if not uniprot_id:
        return None, None
    
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(api_url, timeout=15)
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                pdb_url = data[0].get('pdbUrl')
                pdb_res = requests.get(pdb_url, timeout=20)
                if pdb_res.status_code == 200:
                    return pdb_url, pdb_res.content
    except:
        pass
    return None, None

# --- Streamlit UI 설정 ---
st.set_page_config(page_title="Biobytes Target Analyzer", layout="wide")

st.title("🧬 Biobytes Target Analyzer")
st.markdown("유전자명을 입력하여 **IDG Level(Pharos)**, **임상 현황(Open Targets)**, **구조 데이터(AlphaFold)**를 통합 분석합니다.")

# 입력부
input_text = st.text_input("유전자 기호 입력 (예: ETS2, EGFR)", placeholder="쉼표(,)로 구분하여 여러 개 입력 가능")

if st.button("데이터 분석 및 PDB 찾기"):
    if input_text:
        gene_list = [g.strip().upper() for g in input_text.split(",") if g.strip()]
        
        with st.spinner('데이터를 불러오는 중입니다...'):
            pharos_info = get_pharos_data(gene_list)
            
            if not pharos_info:
                st.error("데이터를 가져오지 못했습니다. 네트워크 상태를 확인하세요.")
            
            for gene in gene_list:
                info = pharos_info.get(gene, {})
                
                # 에러 처리
                if 'error' in info:
                    st.warning(f"⚠️ **{gene}**: {info['error']}")
                    continue

                uniprot_id = info.get('uniprot')
                ot_data = get_opentargets_data(uniprot_id)
                pdb_url, pdb_content = get_alphafold_pdb(uniprot_id)

                # 개별 유전자 리포트 컨테이너
                with st.container():
                    st.write(f"### Target: {gene}")
                    
                    # 상단 요약 지표
                    m_col1, m_col2, m_col3, m_col4 = st.columns(4)
                    with m_col1:
                        st.metric("Full Name", info.get('name', 'N/A'))
                    with m_col2:
                        st.metric("Family", info.get('fam', 'N/A'))
                    with m_col3:
                        st.metric("TDL (개발단계)", info.get('tdl', 'N/A'))
                    with m_col4:
                        disease_count = ot_data.get('associatedDiseases', {}).get('count', 0) if ot_data else 0
                        st.metric("질병 연관성", f"{disease_count} 건")

                    # 상세 섹션 (Expander 사용으로 안정성 확보)
                    with st.expander(f"🔍 {gene} 상세 분석 데이터 (UniProt: {uniprot_id})", expanded=True):
                        
                        col_left, col_right = st.columns(2)
                        
                        with col_left:
                            st.subheader("💊 약물 및 임상 현황")
                            if ot_data and ot_data.get('knownDrugs', {}).get('count', 0) > 0:
                                drugs = ot_data['knownDrugs']['rows']
                                drug_df = pd.DataFrame([
                                    {"Drug Name": d['drug']['name'], "Phase": d['phase'], "Status": d['status']}
                                    for d in drugs[:10]
                                ])
                                # index 오류 방지를 위한 고유 key 부여
                                st.dataframe(drug_df, use_container_width=True, key=f"df_{gene}_{uniprot_id}")
                            else:
                                st.info("알려진 임상 약물 정보가 없습니다.")

                        with col_right:
                            st.subheader("🔬 AlphaFold 구조")
                            if pdb_url and pdb_content:
                                st.success(f"PDB 파일을 찾았습니다.")
                                st.download_button(
                                    label=f"{gene} PDB 다운로드",
                                    data=pdb_content,
                                    file_name=f"AF_{gene}_{uniprot_id}.pdb",
                                    mime="application/octet-stream",
                                    key=f"dl_{gene}_{uniprot_id}" # 유니크한 키 부여
                                )
                                st.caption(f"Source: [AlphaFold EBI]({pdb_url})")
                            else:
                                st.error("AlphaFold PDB를 불러올 수 없습니다.")

                        st.markdown("---")
                        st.subheader("📚 최근 관련 논문")
                        pubs = info.get('publications', [])
                        if pubs:
                            for pub in pubs:
                                date_str = str(pub['date'])[:4] if pub.get('date') else "N/A"
                                st.markdown(f"- **({date_str})** {pub['title']}  \n  *Journal: {pub['journal']}* | [PMID: {pub['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{pub['pmid']}/)")
                        else:
                            st.info("논문 정보가 없습니다.")

                st.divider() # 유전자 간 구분선
    else:
        st.warning("분석할 유전자 기호를 입력해 주세요.")

st.caption("Integrated by Biobytes | Data from Pharos, Open Targets & AlphaFold DB")
