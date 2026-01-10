import streamlit as st
import requests
import pandas as pd

# --- Pharos API 호출 함수 ---
def get_pharos_target_levels(gene_symbols):
    url = 'https://pharos-api.ncats.io/graphql'
    results = []

    query = """
    query getTarget($gene: String!) {
      target(q: { sym: $gene }) {
        sym
        tdl
      }
    }
    """

    for gene in gene_symbols:
        variables = {'gene': gene}
        try:
            response = requests.post(url, json={'query': query, 'variables': variables}, timeout=10)
            if response.status_code == 200:
                data = response.json()
                target_data = data.get('data', {}).get('target')
                if target_data:
                    results.append(target_data)
                else:
                    results.append({'sym': gene, 'tdl': 'Not Found'})
            else:
                results.append({'sym': gene, 'tdl': f'Error ({response.status_code})'})
        except Exception as e:
            results.append({'sym': gene, 'tdl': f'Connection Error'})

    return results

# --- Streamlit UI 설정 ---
st.set_page_config(page_title="Pharos Gene Explorer", layout="centered")

st.title("PharosIDG에서 레벨을 뽑아봅시다.")
st.markdown("유전자 기호를 입력 하면 Level(tclin, tchem 등..)")

# 입력 영역
input_text = st.text_input(
    "유전자 기호 입력 (쉼표로 구분)", 
    placeholder="예: EGFR, TP53, BRCA1",
    help="여러 개를 조회할 경우 쉼표(,)로 구분하여 입력하세요."
)

if st.button("데이터 조회하기"):
    if input_text:
        # 1. 입력 문자열 처리 (쉼표 분리 -> 공백 제거 -> 대문자 변환)
        gene_list = [g.strip().upper() for g in input_text.split(",") if g.strip()]
        
        with st.spinner('데이터를 가져오는 중입니다...'):
            # 2. API 호출
            data = get_pharos_target_levels(gene_list)
            
            # 3. 데이터프레임 변환 및 표시
            df = pd.DataFrame(data)
            
            # 컬럼명 변경
            df.columns = ["Gene Symbol", "IDG Level (TDL)"]
            
            st.divider() # 구분선
            
            # 결과 표 출력
            st.subheader("조회 결과")
            st.dataframe(df, use_container_width=True, hide_index=True)
            
            # 다운로드 버튼 추가 (CSV)
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="결과를 CSV로 저장",
                data=csv,
                file_name="pharos_results.csv",
                mime="text/csv",
            )
    else:
        st.warning("조회할 유전자 기호를 입력해주세요.")

# 하단 정보
st.caption("Data source: Pharos (NCATS)")