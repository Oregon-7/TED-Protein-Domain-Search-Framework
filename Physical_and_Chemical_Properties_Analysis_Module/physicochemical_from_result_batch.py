import re
from pathlib import Path

from physicochemical_properties_batch import run_physicochemical_pipeline


CATH_PATTERN = re.compile(r"^\d+\.\d+\.\d+(?:\.\d+)?$")
SKIP_PREFIXES = (
    "domain_sequences_",
    "domain_properties_",
    "domain_failed_records_",
    "~$",
)


def is_source_input_file(file_path: Path) -> bool:
    name = file_path.name
    if any(name.startswith(prefix) for prefix in SKIP_PREFIXES):
        return False

    stem = file_path.stem
    if "_" not in stem:
        return False

    cath_candidate = stem.split("_", 1)[0]
    return bool(CATH_PATTERN.match(cath_candidate))


def extract_cath_from_filename(file_path: Path) -> str:
    return file_path.stem.split("_", 1)[0]


def collect_source_files(result_root: Path) -> list[Path]:
    files = [p for p in result_root.rglob("*.xlsx") if is_source_input_file(p)]
    return sorted(files)


def main() -> None:
    root_input = input("请输入Result根目录路径（例如 C:\\Users\\C\\PythonBatch\\workflow_260420\\Result）: ").strip().strip('"')
    result_root = Path(root_input)

    if not result_root.exists() or not result_root.is_dir():
        print(f"错误：目录不存在或不是文件夹: {result_root}")
        return

    source_files = collect_source_files(result_root)
    if not source_files:
        print("未找到可处理的来源xlsx文件。")
        return

    print(f"检测到 {len(source_files)} 个待处理文件。")

    success = 0
    failed = 0

    for idx, file_path in enumerate(source_files, 1):
        cath_id = extract_cath_from_filename(file_path)
        output_dir = file_path.parent

        print("-" * 60)
        print(f"[{idx}/{len(source_files)}] 文件: {file_path}")
        print(f"  解析CATH: {cath_id}")
        print(f"  输出目录: {output_dir}")

        try:
            result = run_physicochemical_pipeline(
                excel_path=str(file_path),
                target_cath=cath_id,
                output_dir=output_dir,
            )
            if result is None:
                failed += 1
                print("  结果: 失败")
            else:
                success += 1
                print("  结果: 成功")
        except Exception as exc:
            failed += 1
            print(f"  结果: 异常失败 ({exc})")

    print("=" * 60)
    print("批量分析完成")
    print(f"成功: {success}")
    print(f"失败: {failed}")


if __name__ == "__main__":
    main()
